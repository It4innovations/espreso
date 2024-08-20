
#ifdef HAVE_CUDA
#ifndef USE_CUSPARSE_LEGACY

#include "gpu/gpu_spblas.h"
#include "w.cuda.gpu_management.h"
#include "basis/utilities/utils.h"

#include <cusparse.h>
#include <complex>



inline void _check(cusparseStatus_t status, const char *file, int line)
{
    if (status != CUSPARSE_STATUS_SUCCESS) {
        espreso::eslog::error("CUSPARSE Error %d %s: %s. In file '%s' on line %d\n", status, cusparseGetErrorName(status), cusparseGetErrorString(status), file, line);
    }
}



#include "w.cuda.gpu_spblas_common.h"



namespace espreso {
namespace gpu {
namespace spblas {

    spblas_wrapper_impl get_implementation()
    {
        return spblas_wrapper_impl::CUSPARSE_MODERN;
    }

    namespace
    {
        template<typename I>
        static cusparseIndexType_t _sparse_index_type()
        {
            if constexpr(std::is_same_v<I, int32_t>) return CUSPARSE_INDEX_32I;
            if constexpr(std::is_same_v<I, int64_t>) return CUSPARSE_INDEX_64I;
        }

        template<typename T>
        static cudaDataType_t _sparse_data_type()
        {
            if constexpr(std::is_same_v<T, float>)  return CUDA_R_32F;
            if constexpr(std::is_same_v<T, double>) return CUDA_R_64F;
            if constexpr(std::is_same_v<T, std::complex<float>>)  return CUDA_C_32F;
            if constexpr(std::is_same_v<T, std::complex<double>>) return CUDA_C_64F;
        }

        template<typename T>
        static cusparseOperation_t _char_to_operation(char c)
        {
            if(utils::is_real<T>()) {
                switch(c) {
                    case 'N': return CUSPARSE_OPERATION_NON_TRANSPOSE;
                    case 'T': return CUSPARSE_OPERATION_TRANSPOSE;
                    case 'H': return CUSPARSE_OPERATION_TRANSPOSE;
                    default: eslog::error("invalid operation '%c'\n", c);
                }
            }
            if(utils::is_complex<T>()) {
                switch(c) {
                    case 'N': return CUSPARSE_OPERATION_NON_TRANSPOSE;
                    case 'T': return CUSPARSE_OPERATION_TRANSPOSE;
                    case 'H': return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
                    default: eslog::error("invalid operation '%c'\n", c);
                }
            }
            eslog::error("invalid type\n");
        }
    }

    struct _handle
    {
        cusparseHandle_t h;
        cudaStream_t get_stream()
        {
            cudaStream_t stream;
            CHECK(cusparseGetStream(h, &stream));
            return stream;
        }
    };

    struct _descr_matrix_csr
    {
        cusparseSpMatDescr_t d;
        cusparseSpMatDescr_t d_for_sp2dn; // no attributes set
        void * vals_ptr;
    };

    struct _descr_matrix_dense
    {
        cusparseDnMatDescr_t d;
        cusparseDnMatDescr_t d_complementary;
        _descr_matrix_dense get_complementary()
        {
            _descr_matrix_dense ret;
            ret.d = d_complementary;
            ret.d_complementary = d;
            return ret;
        }
    };

    struct _descr_vector_dense
    {
        cusparseDnVecDescr_t d;
    };

    struct _descr_sparse_trsv
    {
        cusparseSpSVDescr_t d;
    };

    struct _descr_sparse_trsm
    {
        cusparseSpSMDescr_t d;
    };

    struct _descr_sparse_mv
    {
#if CUDART_VERSION >= 12040 && CUDART_VERSION < 12070 // see the mv function for info
        cusparseSpMatDescr_t copy_of_matrix_descr;
        bool was_matrix_descr_initialized = false;
#endif
    };

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        CHECK(cusparseCreate(&h->h));
        CHECK(cusparseSetStream(h->h, q->stream));
    }

    void handle_destroy(handle & h)
    {
        if(h.get() == nullptr) return;

        CHECK(cusparseDestroy(h->h));
        h.reset();
    }

    template<typename T, typename I>
    void descr_matrix_csr_create(handle & /*h*/, descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill)
    {
        descr = std::make_shared<_descr_matrix_csr>();
        void * dummyptr = reinterpret_cast<void*>(8);
        CHECK(cusparseCreateCsr(&descr->d_for_sp2dn, nrows, ncols, nnz, dummyptr, dummyptr, dummyptr, _sparse_index_type<I>(), _sparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, _sparse_data_type<T>()));
        CHECK(cusparseCreateCsr(&descr->d,           nrows, ncols, nnz, dummyptr, dummyptr, dummyptr, _sparse_index_type<I>(), _sparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, _sparse_data_type<T>()));
        auto upper = CUSPARSE_FILL_MODE_UPPER;
        auto lower = CUSPARSE_FILL_MODE_LOWER;
        auto nonunit = CUSPARSE_DIAG_TYPE_NON_UNIT;
        if(fill == 'L') CHECK(cusparseSpMatSetAttribute(descr->d, CUSPARSE_SPMAT_FILL_MODE, &lower, sizeof(lower)));
        if(fill == 'U') CHECK(cusparseSpMatSetAttribute(descr->d, CUSPARSE_SPMAT_FILL_MODE, &upper, sizeof(upper)));
        if(fill != 'N') CHECK(cusparseSpMatSetAttribute(descr->d, CUSPARSE_SPMAT_DIAG_TYPE, &nonunit, sizeof(nonunit)));
    }

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(handle & /*h*/, descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(cusparseCsrSetPointers(descr->d,           matrix.rows, matrix.cols, matrix.vals));
        CHECK(cusparseCsrSetPointers(descr->d_for_sp2dn, matrix.rows, matrix.cols, matrix.vals));
        descr->vals_ptr = matrix.vals;
    }

    void descr_matrix_csr_destroy(handle & /*h*/, descr_matrix_csr & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(cusparseDestroySpMat(descr->d));
        CHECK(cusparseDestroySpMat(descr->d_for_sp2dn));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_matrix_dense_create(handle & /*h*/, descr_matrix_dense & descr, I nrows, I ncols, I ld, char order)
    {
        descr = std::make_shared<_descr_matrix_dense>();
        void * dummyptr = reinterpret_cast<void*>(8);
        if(order == 'R') CHECK(cusparseCreateDnMat(&descr->d, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_ROW));
        if(order == 'C') CHECK(cusparseCreateDnMat(&descr->d, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_COL));
        if(order == 'R') CHECK(cusparseCreateDnMat(&descr->d_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_COL));
        if(order == 'C') CHECK(cusparseCreateDnMat(&descr->d_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_ROW));
    }

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(handle & /*h*/, descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(cusparseDnMatSetValues(descr->d, matrix.vals));
        CHECK(cusparseDnMatSetValues(descr->d_complementary, matrix.vals));
    }

    void descr_matrix_dense_destroy(handle & /*h*/, descr_matrix_dense & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(cusparseDestroyDnMat(descr->d));
        CHECK(cusparseDestroyDnMat(descr->d_complementary));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_vector_dense_create(handle & /*h*/, descr_vector_dense & descr, I size)
    {
        descr = std::make_shared<_descr_vector_dense>();
        void * dummyptr = reinterpret_cast<void*>(8);
        CHECK(cusparseCreateDnVec(&descr->d, size, dummyptr, _sparse_data_type<T>()));
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & /*h*/, descr_vector_dense & descr, Vector_Dense<T,I,A> & vector)
    {
        static_assert(A::is_data_device_accessible, "vector data must be device accessible");
        CHECK(cusparseDnVecSetValues(descr->d, vector.vals));
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & /*h*/, descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx)
    {
        static_assert(A::is_data_device_accessible, "vector data must be device accessible");
        CHECK(cusparseDnVecSetValues(descr->d, matrix.vals + colidx * matrix.get_ld()));
    }

    void descr_vector_dense_destroy(handle & /*h*/, descr_vector_dense & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(cusparseDestroyDnVec(descr->d));
        descr.reset();
    }

    void descr_sparse_trsv_create(handle & /*h*/, descr_sparse_trsv & descr)
    {
        descr = std::make_shared<_descr_sparse_trsv>();
        CHECK(cusparseSpSV_createDescr(&descr->d));
    }

    void descr_sparse_trsv_destroy(handle & /*h*/, descr_sparse_trsv & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(cusparseSpSV_destroyDescr(descr->d));
        descr.reset();
    }

    void descr_sparse_trsm_create(handle & /*h*/, descr_sparse_trsm & descr)
    {
        descr = std::make_shared<_descr_sparse_trsm>();
        CHECK(cusparseSpSM_createDescr(&descr->d));
    }

    void descr_sparse_trsm_destroy(handle & /*h*/, descr_sparse_trsm & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(cusparseSpSM_destroyDescr(descr->d));
        descr.reset();
    }

    void descr_sparse_mv_create(handle & /*h*/, descr_sparse_mv & descr)
    {
#if CUDART_VERSION >= 12040 && CUDART_VERSION < 12070 // see the mv function for info
        descr = std::make_shared<_descr_sparse_mv>();
#endif
    }

    void descr_sparse_mv_destroy(handle & /*h*/, descr_sparse_mv & descr)
    {
#if CUDART_VERSION >= 12040 && CUDART_VERSION < 12070 // see the mv function for info
        if(descr.get() == nullptr) return;

        if(descr->was_matrix_descr_initialized)
        {
            CHECK(cusparseDestroySpMat(descr->copy_of_matrix_descr));
        }
        descr.reset();
#endif
    }

    template<typename T, typename I>
    void transpose(handle & h, descr_matrix_csr & output, descr_matrix_csr & input, bool conjugate, size_t & buffersize, void * buffer, char stage)
    {
        int64_t out_nrows, out_ncols, out_nnz, in_nrows, in_ncols, in_nnz;
        void *out_rowptrs, *out_colidxs, *out_vals, *in_rowptrs, *in_colidxs, *in_vals;
        cusparseIndexType_t rowoffsettype, colindtype;
        cusparseIndexBase_t idxbase;
        cudaDataType type;
        CHECK(cusparseCsrGet(output->d, &out_nrows, &out_ncols, &out_nnz, &out_rowptrs, &out_colidxs, &out_vals, &rowoffsettype, &colindtype, &idxbase, &type));
        CHECK(cusparseCsrGet(input->d,  &in_nrows,  &in_ncols,  &in_nnz,  &in_rowptrs,  &in_colidxs,  &in_vals,  &rowoffsettype, &colindtype, &idxbase, &type));
        cudaStream_t stream = h->get_stream();
        if(stage == 'B') my_csr_transpose_buffersize<I>(stream, in_nrows, in_ncols, in_nnz, buffersize);
        if(stage == 'P') my_csr_transpose_preprocess<I>(stream, in_nrows, in_ncols, in_nnz, (I*)in_rowptrs, (I*)in_colidxs, (I*)out_rowptrs, (I*)out_colidxs, buffersize, buffer);
        if(stage == 'C') my_csr_transpose_compute<T,I>(stream, in_nnz, (T*)in_vals, (T*)out_vals, conjugate, buffer);
    }

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage)
    {
        if(transpose == 'N') {
            if(stage == 'B') CHECK(cusparseSparseToDense_bufferSize(h->h, sparse->d_for_sp2dn, dense->d, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, &buffersize));
            if(stage == 'C') CHECK(cusparseSparseToDense           (h->h, sparse->d_for_sp2dn, dense->d, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, buffer));
        }
        else if(transpose == 'T') {
            descr_matrix_dense descr_dense_complementary = std::make_shared<_descr_matrix_dense>(dense->get_complementary());
            sparse_to_dense<T,I>(h, 'N', sparse, descr_dense_complementary, buffersize, buffer, stage);
            return;
        }
        else {
            eslog::error("transpose '%c' not supported\n", transpose);
        }
    }

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage)
    {
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpSV_bufferSize  (h->h, _char_to_operation<T>(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, &buffersize));
#if CUDART_VERSION >= 12020 // cusparseSpSV_updateMatrix available since CUDA/12.1.1, but no way to check for the .1 update
        if(stage == 'P') CHECK(cusparseSpSV_analysis    (h->h, _char_to_operation<T>(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, buffer));
        if(stage == 'U') CHECK(cusparseSpSV_updateMatrix(h->h, descr_trsv->d, matrix->vals_ptr, CUSPARSE_SPSV_UPDATE_GENERAL));
#else
        // if(stage == 'P') CHECK(cusparseSpSV_analysis    (h->h, _char_to_operation<T>(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, buffer)); // no need to do preprocess when update must happen anyway
        if(stage == 'U') CHECK(cusparseSpSV_analysis    (h->h, _char_to_operation<T>(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, buffer));
#endif
        if(stage == 'C') CHECK(cusparseSpSV_solve       (h->h, _char_to_operation<T>(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d));
    }

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
    {
        if(transpose_sol == 'N') {
            T one = 1.0;
            if(stage == 'B') CHECK(cusparseSpSM_bufferSize  (h->h, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, &buffersize));
#if CUDART_VERSION >= 12040
            if(stage == 'P') CHECK(cusparseSpSM_analysis    (h->h, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, buffer));
            if(stage == 'U') CHECK(cusparseSpSM_updateMatrix(h->h, descr_trsm->d, matrix->vals_ptr, CUSPARSE_SPSM_UPDATE_GENERAL));
#else
            // if(stage == 'P') CHECK(cusparseSpSM_analysis    (h->h, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, buffer)); // it does not make sense to do Preprocessing. Update will be called anyway, which has to completely redo the analysis
            if(stage == 'U') CHECK(cusparseSpSM_analysis    (h->h, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, buffer));
#endif
            if(stage == 'C') CHECK(cusparseSpSM_solve       (h->h, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d));
        }
        else if(transpose_sol == 'T') {
            descr_matrix_dense descr_sol_compl = std::make_shared<_descr_matrix_dense>(sol->get_complementary());
            trsm<T,I>(h, transpose_mat, transpose_rhs, 'N', matrix, rhs, sol, descr_trsm, buffersize, buffer, stage);
        }
        else {
            eslog::error("transpose_sol '%c' not supported\n", transpose_sol);
        }
    }

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, descr_sparse_mv & descr_mv, size_t & buffersize, void * buffer, char stage)
    {
        cusparseSpMatDescr_t * Ad = &A->d;
#if CUDART_VERSION >= 12040 && CUDART_VERSION < 12070 // hopefully fixed in 12.7, did not check
        {
            // https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html#cusparse-release-12-6
            // cannot perform different operations with the same handle and different buffer
            // either use the same workspace buffer
            // or use separate matrix handle for each spmv operation (this solution)
            int64_t nrows, ncols, nnz;
            void *rowptrs, *colidxs, *vals;
            cusparseIndexType_t ro_type, ci_type;
            cusparseIndexBase_t idxbase;
            cudaDataType value_type;
            CHECK(cusparseCsrGet(A->d, &nrows, &ncols, &nnz, &rowptrs, &colidxs, &vals, &ro_type, &ci_type, &idxbase, &value_type));

            if(!descr_mv->was_matrix_descr_initialized)
            {
                CHECK(cusparseCreateCsr(&descr_mv->copy_of_matrix_descr, nrows, ncols, nnz, rowptrs, colidxs, vals, ro_type, ci_type, idxbase, value_type));
                descr_mv->was_matrix_descr_initialized = true;
            }
            CHECK(cusparseCsrSetPointers(descr_mv->copy_of_matrix_descr, rowptrs, colidxs, vals));
            Ad = &descr_mv->copy_of_matrix_descr;
        }
#endif

        T one = 1.0;
        T zero = 0.0;
        if(stage == 'B') CHECK(cusparseSpMV_bufferSize(h->h, _char_to_operation<T>(transpose), &one, *Ad, x->d, &zero, y->d, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, &buffersize));
#if CUDART_VERSION >= 12040
        if(stage == 'P') CHECK(cusparseSpMV_preprocess(h->h, _char_to_operation<T>(transpose), &one, *Ad, x->d, &zero, y->d, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, buffer));
#else
        // if(stage == 'P') ; // no preprocess function exists in CUDA < 12.4
#endif
        if(stage == 'C') CHECK(cusparseSpMV           (h->h, _char_to_operation<T>(transpose), &one, *Ad, x->d, &zero, y->d, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, buffer));
    }

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpMM_bufferSize(h->h, _char_to_operation<T>(transpose_A), _char_to_operation<T>(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, &buffersize));
        if(stage == 'P') CHECK(cusparseSpMM_preprocess(h->h, _char_to_operation<T>(transpose_A), _char_to_operation<T>(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
        if(stage == 'C') CHECK(cusparseSpMM           (h->h, _char_to_operation<T>(transpose_A), _char_to_operation<T>(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
    }

}
}
}

#include "gpu/gpu_spblas.inst.hpp"

#endif
#endif
