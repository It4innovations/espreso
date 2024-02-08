
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
        char str[1000];
        snprintf(str, sizeof(str), "CUSPARSE Error %d %s: %s. In file '%s' on line %d\n", status, cusparseGetErrorName(status), cusparseGetErrorString(status), file, line);
        espreso::eslog::error(str);
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
        if(transpose == 'T')
        {
            descr_matrix_dense descr_dense_complementary = dense.get_complementary();
            sparse_to_dense<T,I>(h, 'N', sparse, descr_dense_complementary, buffersize, buffer, stage);
            return;
        }
        if(stage == 'B') CHECK(cusparseSparseToDense_bufferSize(h->h, sparse->d, dense->d, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, &buffersize));
        if(stage == 'C') CHECK(cusparseSparseToDense           (h->h, sparse->d, dense->d, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, buffer));
    }

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage)
    {
        cusparseOperation_t op = (transpose == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpSV_bufferSize  (h->h, op, &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, &buffersize));
        if(stage == 'P') CHECK(cusparseSpSV_analysis    (h->h, op, &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d, buffer));
        if(stage == 'U') CHECK(cusparseSpSV_updateMatrix(h->h, descr_trsv->d, matrix->vals_ptr, CUSPARSE_SPSV_UPDATE_GENERAL));
        if(stage == 'C') CHECK(cusparseSpSV_solve       (h->h, op, &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSV_ALG_DEFAULT, descr_trsv->d));
    }

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
    {
        cusparseOperation_t op_mat = (transpose_mat == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        cusparseOperation_t op_rhs = (transpose_rhs == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpSM_bufferSize(h->h, op_mat, op_rhs, &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, &buffersize));
        // if(stage == 'P') CHECK(cusparseSpSM_analysis  (h->h, op_mat, op_rhs, &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, buffer)); // it does not make sense to do Preprocessing. Update will be called anyway, which has to completely redo the analysis
        if(stage == 'U') CHECK(cusparseSpSM_analysis  (h->h, op_mat, op_rhs, &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d, buffer));
        if(stage == 'C') CHECK(cusparseSpSM_solve     (h->h, op_mat, op_rhs, &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), CUSPARSE_SPSM_ALG_DEFAULT, descr_trsm->d));
    }

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage)
    {
        cusparseOperation_t op = (transpose == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        T one = 1.0;
        T zero = 0.0;
        if(stage == 'B') CHECK(cusparseSpMV_bufferSize(h->h, op, &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, &buffersize));
        if(stage == 'C') CHECK(cusparseSpMV           (h->h, op, &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, buffer));
    }

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        cusparseOperation_t op_A = (transpose_A == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        cusparseOperation_t op_B = (transpose_B == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpMM_bufferSize(h->h, op_A, op_B, &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, &buffersize));
        if(stage == 'P') CHECK(cusparseSpMM_preprocess(h->h, op_A, op_B, &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
        if(stage == 'C') CHECK(cusparseSpMM           (h->h, op_A, op_B, &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
    }



    #define INSTANTIATE(T,I) \
    template void descr_matrix_csr_create<T,I>(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char symmetry); \
    template void descr_matrix_dense_create<T,I>(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order); \
    template void descr_vector_dense_create<T,I>(descr_vector_dense & descr, I size); \
    template void sparse_to_dense<T,I>(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage); \
    template void trsv<T,I>(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage); \
    template void trsm<T,I>(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage); \
    template void mv<T,I>(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage); \
    template void mm<T,I>(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage);
        // INSTANTIATE(float,                int32_t)
        INSTANTIATE(double,               int32_t)
        // INSTANTIATE(std::complex<float >, int32_t)
        // INSTANTIATE(std::complex<double>, int32_t)
        // INSTANTIATE(float,                int64_t)
        // INSTANTIATE(double,               int64_t)
        // INSTANTIATE(std::complex<float >, int64_t)
        // INSTANTIATE(std::complex<double>, int64_t)
    #undef INSTANTIATE

    #define INSTANTIATE(T,I,Adevice) \
    template void descr_matrix_csr_link_data<T,I,Adevice>(descr_matrix_csr & descr, Matrix_CSR<T,I,Adevice> & matrix); \
    template void descr_matrix_dense_link_data<T,I,Adevice>(descr_matrix_dense & descr, Matrix_Dense<T,I,Adevice> & matrix); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Vector_Dense<T,I,Adevice> & vector); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Matrix_Dense<T,I,Adevice> & matrix, I colidx = 0); \
        // INSTANTIATE(float,                int32_t, mgm::Ad)
        INSTANTIATE(double,               int32_t, mgm::Ad)
        // INSTANTIATE(std::complex<float >, int32_t, mgm::Ad)
        // INSTANTIATE(std::complex<double>, int32_t, mgm::Ad)
        // INSTANTIATE(float,                int64_t, mgm::Ad)
        // INSTANTIATE(double,               int64_t, mgm::Ad)
        // INSTANTIATE(std::complex<float >, int64_t, mgm::Ad)
        // INSTANTIATE(std::complex<double>, int64_t, mgm::Ad)
        // INSTANTIATE(float,                int32_t, cbmba_d)
        INSTANTIATE(double,               int32_t, cbmba_d)
        // INSTANTIATE(std::complex<float >, int32_t, cbmba_d)
        // INSTANTIATE(std::complex<double>, int32_t, cbmba_d)
        // INSTANTIATE(float,                int64_t, cbmba_d)
        // INSTANTIATE(double,               int64_t, cbmba_d)
        // INSTANTIATE(std::complex<float >, int64_t, cbmba_d)
        // INSTANTIATE(std::complex<double>, int64_t, cbmba_d)
    #undef INSTANTIATE

}
}
}

#endif
#endif
