
#ifdef HAVE_CUDA
#ifndef USE_CUSPARSE_LEGACY

#include "gpu/gpu_spblas.h"
#include "w.cuda.gpu_management.h"

#include <cusparse.h>
#include <complex>



inline void _check(cusparseStatus_t status, const char *file, int line)
{
    if (status != CUSPARSE_STATUS_SUCCESS)
    {
        espreso::eslog::error("CUSPARSE Error %d %s: %s. In file '%s' on line %d\n", status, cusparseGetErrorName(status), cusparseGetErrorString(status), file, line);
    }
}



namespace espreso {
namespace gpu {
namespace spblas {

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

        static cusparseOperation_t _char_to_operation(char c)
        {
            switch(c)
            {
                case 'N': return CUSPARSE_OPERATION_TRANSPOSE;
                case 'T': return CUSPARSE_OPERATION_TRANSPOSE;
                case 'H': return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
                default: eslog::error("invalid operation '%c'\n", c);
            }
        }
    }

    struct _handle
    {
        cusparseHandle_t h;
    };

    struct _descr_matrix_csr
    {
        cusparseSpMatDescr_t d;
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

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        CHECK(cusparseCreate(&h->h));
        CHECK(cusparseSetStream(h->h, q->stream));
    }

    void handle_destroy(handle & h)
    {
        CHECK(cusparseDestroy(h->h));
        h.reset();
    }

    template<typename T, typename I>
    void descr_matrix_csr_create(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char symmetry)
    {
        descr = std::make_shared<_descr_matrix_csr>();
        void * dummyptr = reinterpret_cast<void*>(1);
        CHECK(cusparseCreateCsr(&descr->d, nrows, ncols, nnz, dummyptr, dummyptr, dummyptr, _sparse_index_type<I>(), _sparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, _sparse_data_type<T>()));
        auto upper = CUSPARSE_FILL_MODE_UPPER;
        auto lower = CUSPARSE_FILL_MODE_LOWER;
        auto nonunit = CUSPARSE_DIAG_TYPE_NON_UNIT;
        if(symmetry == 'L') CHECK(cusparseSpMatSetAttribute(descr->d, CUSPARSE_SPMAT_FILL_MODE, &lower, sizeof(lower)));
        if(symmetry == 'U') CHECK(cusparseSpMatSetAttribute(descr->d, CUSPARSE_SPMAT_FILL_MODE, &upper, sizeof(upper)));
        if(symmetry != 'N') CHECK(cusparseSpMatSetAttribute(descr->d, CUSPARSE_SPMAT_DIAG_TYPE, &nonunit, sizeof(nonunit)));
    }

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(cusparseCsrSetPointers(descr->d, matrix.rows, matrix.cols, matrix.vals));
        descr->vals_ptr = matrix.vals;
    }

    void descr_matrix_csr_destroy(descr_matrix_csr & descr)
    {
        CHECK(cusparseDestroySpMat(descr->d));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_matrix_dense_create(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order)
    {
        descr = std::make_shared<_descr_matrix_dense>();
        void * dummyptr = reinterpret_cast<void*>(1);
        if(order == 'R') CHECK(cusparseCreateDnMat(&descr->d, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_ROW));
        if(order == 'C') CHECK(cusparseCreateDnMat(&descr->d, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_COL));
        if(order == 'R') CHECK(cusparseCreateDnMat(&descr->d_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_COL));
        if(order == 'C') CHECK(cusparseCreateDnMat(&descr->d_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_ROW));
    }

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(cusparseDnMatSetValues(descr->d, matrix.vals));
        CHECK(cusparseDnMatSetValues(descr->d_complementary, matrix.vals));
    }

    void descr_matrix_dense_destroy(descr_matrix_dense & descr)
    {
        CHECK(cusparseDestroyDnMat(descr->d));
        CHECK(cusparseDestroyDnMat(descr->d_complementary));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_vector_dense_create(descr_vector_dense & descr, I size)
    {
        descr = std::make_shared<_descr_vector_dense>();
        void * dummyptr = reinterpret_cast<void*>(1);
        CHECK(cusparseCreateDnVec(&descr->d, size, dummyptr, _sparse_data_type<T>()));
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<T,I,A> & vector)
    {
        static_assert(A::is_data_device_accessible, "vector data must be device accessible");
        CHECK(cusparseDnVecSetValues(descr->d, vector.vals));
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx)
    {
        static_assert(A::is_data_device_accessible, "vector data must be device accessible");
        CHECK(cusparseDnVecSetValues(descr->d, matrix.vals + colidx * matrix.get_ld()));
    }

    void descr_vector_dense_destroy(descr_vector_dense & descr)
    {
        CHECK(cusparseDestroyDnVec(descr->d));
        descr.reset();
    }

    void descr_sparse_trsv_create(descr_sparse_trsv & descr)
    {
        descr = std::make_shared<_descr_sparse_trsv>();
        CHECK(cusparseSpSV_createDescr(&descr->d));
    }

    void descr_sparse_trsv_destroy(descr_sparse_trsv & descr)
    {
        CHECK(cusparseSpSV_destroyDescr(descr->d));
        descr.reset();
    }

    void descr_sparse_trsm_create(descr_sparse_trsm & descr)
    {
        descr = std::make_shared<_descr_sparse_trsm>();
        CHECK(cusparseSpSM_createDescr(&descr->d));
    }

    void descr_sparse_trsm_destroy(descr_sparse_trsm & descr)
    {
        CHECK(cusparseSpSM_destroyDescr(descr->d));
        descr.reset();
    }

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage)
    {
        if(transpose == 'N')
        {
            if(stage == 'B') CHECK(cusparseSparseToDense_bufferSize(h->h, sparse->d, dense->d, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, &buffersize));
            if(stage == 'C') CHECK(cusparseSparseToDense           (h->h, sparse->d, dense->d, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, buffer));
        }
        else if(transpose == 'T')
        {
            descr_matrix_dense descr_dense_complementary = std::make_shared<_descr_matrix_dense>(dense->get_complementary());
            sparse_to_dense<T,I>(h, 'N', sparse, descr_dense_complementary, buffersize, buffer, stage);
            return;
        }
        else
        {
            eslog::error("transpose '%c' not supported\n", transpose);
        }
    }

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage)
    {
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpSV_bufferSize  (h->h, _char_to_operation(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, &buffersize));
        if(stage == 'P') CHECK(cusparseSpSV_analysis    (h->h, _char_to_operation(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, buffer));
        if(stage == 'U') CHECK(cusparseSpSV_updateMatrix(h->h, descr_trsv->d, matrix->vals_ptr, CUSPARSE_SPSV_UPDATE_GENERAL));
        if(stage == 'C') CHECK(cusparseSpSV_solve       (h->h, _char_to_operation(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d));
    }

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
    {
        if(transpose_sol == 'N')
        {
            T one = 1.0;
            if(stage == 'B') CHECK(cusparseSpSM_bufferSize(h->h, _char_to_operation(transpose_mat), _char_to_operation(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, &buffersize));
            // if(stage == 'P') CHECK(cusparseSpSM_analysis  (h->h, _char_to_operation(transpose_mat), _char_to_operation(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, buffer)); // it does not make sense to do Preprocessing. Update will be called anyway, which has to completely redo the analysis
            if(stage == 'U') CHECK(cusparseSpSM_analysis  (h->h, _char_to_operation(transpose_mat), _char_to_operation(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, buffer));
            if(stage == 'C') CHECK(cusparseSpSM_solve     (h->h, _char_to_operation(transpose_mat), _char_to_operation(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d));
        }
        else if(transpose_sol == 'T')
        {
            descr_matrix_dense descr_sol_compl = std::make_shared<_descr_matrix_dense>(sol->get_complementary());
            trsm(h, transpose_mat, transpose_rhs, 'N', matrix, rhs, sol, descr_trsm, buffersize, buffer, stage);
        }
        else
        {
            eslog::error("transpose_sol '%c' not supported\n", transpose);
        }
    }

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage)
    {
        T one = 1.0;
        T zero = 0.0;
        if(stage == 'B') CHECK(cusparseSpMV_bufferSize(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, &buffersize));
        if(stage == 'C') CHECK(cusparseSpMV           (h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, buffer));
    }

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpMM_bufferSize(h->h, _char_to_operation(transpose_A), _char_to_operation(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, &buffersize));
        if(stage == 'P') CHECK(cusparseSpMM_preprocess(h->h, _char_to_operation(transpose_A), _char_to_operation(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
        if(stage == 'C') CHECK(cusparseSpMM           (h->h, _char_to_operation(transpose_A), _char_to_operation(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
    }



    #define INSTANTIATE_T_I_ADEVICE(T,I,Adevice) \
    template void descr_matrix_csr_link_data<T,I,Adevice>(descr_matrix_csr & descr, Matrix_CSR<T,I,Adevice> & matrix); \
    template void descr_matrix_dense_link_data<T,I,Adevice>(descr_matrix_dense & descr, Matrix_Dense<T,I,Adevice> & matrix); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Vector_Dense<T,I,Adevice> & vector); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Matrix_Dense<T,I,Adevice> & matrix, I colidx = 0);

        #define INSTANTIATE_T_I(T,I) \
        template void descr_matrix_csr_create<T,I>(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char symmetry); \
        template void descr_matrix_dense_create<T,I>(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order); \
        template void descr_vector_dense_create<T,I>(descr_vector_dense & descr, I size); \
        template void sparse_to_dense<T,I>(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage); \
        template void trsv<T,I>(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage); \
        template void trsm<T,I>(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage); \
        template void mv<T,I>(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage); \
        template void mm<T,I>(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage); \
        INSTANTIATE_T_I_ADEVICE(T, I, mgm::Ad) \
        INSTANTIATE_T_I_ADEVICE(T, I, cbmba_d)

            #define INSTANTIATE_T(T) \
            INSTANTIATE_T_I(T, int32_t) \
            /* INSTANTIATE_T_I(T, int64_t) */

                // INSTANTIATE_T(float)
                INSTANTIATE_T(double)
                // INSTANTIATE_T(std::complex<float>)
                // INSTANTIATE_T(std::complex<double>)

            #undef INSTANTIATE_T
        #undef INSTANTIATE_T_I
    #undef INSTANTIATE_T_I_ADEVICE

}
}
}

#endif
#endif
