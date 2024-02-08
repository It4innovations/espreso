
#ifdef HAVE_CUDA
#ifdef USE_CUSPARSE_LEGACY

#include "gpu/gpu_spblas.h"
#include "w.cuda.gpu_management.h"

#include <cusparse.h>
#include <complex>

#if defined(__GNUC__) && !defined(__clang__)
#define MY_COMPILER_GCC
#elif defined(__clang__)
#define MY_COMPILER_CLANG
#endif

// The legacy cusparse API is deprecated. I know, no need to remind me.
#ifdef MY_COMPILER_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#ifdef MY_COMPILER_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif



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
        template<typename T> struct cpp_to_cuda_type { using type = T; };
        template<> struct cpp_to_cuda_type<std::complex<float>> { using type = cuComplex; };
        template<> struct cpp_to_cuda_type<std::complex<double>> { using type = cuDoubleComplex; };
        template<typename T> using cpp_to_cuda_type_t = typename cpp_to_cuda_type<T>::type;

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

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsv_buffersize(cusparseHandle_t handle, cusparseOperation_t transA, I m, I nnz, const cusparseMatDescr_t descrA, void * csrValA, const void * csrRowPtrA, const void * csrColIndA, csrsv2Info_t info, int * pBufferSizeInBytes)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, pBufferSizeInBytes);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, pBufferSizeInBytes);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, pBufferSizeInBytes);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, pBufferSizeInBytes);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsv_analysis(cusparseHandle_t handle, cusparseOperation_t transA, I m, I nnz, const cusparseMatDescr_t descrA, void * csrValA, const void * csrRowPtrA, const void * csrColIndA, csrsv2Info_t info, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, policy, pBuffer);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsv_solve(cusparseHandle_t handle, cusparseOperation_t transA, I m, I nnz, void * alpha, const cusparseMatDescr_t descrA, void * csrValA, const void * csrRowPtrA, const void * csrColIndA, csrsv2Info_t info, const void * x, void * y, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsv2_solve(handle, transA, m, nnz, (U*)alpha, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, (U*)x, (U*)y, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsv2_solve(handle, transA, m, nnz, (U*)alpha, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, (U*)x, (U*)y, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsv2_solve(handle, transA, m, nnz, (U*)alpha, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, (U*)x, (U*)y, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsv2_solve(handle, transA, m, nnz, (U*)alpha, descrA, (U*)csrValA, (I*)csrRowPtrA, (I*)csrColIndA, info, (U*)x, (U*)y, policy, pBuffer);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsm_buffersize(cusparseHandle_t handle, int algo, cusparseOperation_t transA, cusparseOperation_t transB, int m, int nrhs, int nnz, const void * alpha, const cusparseMatDescr_t descrA, const void * csrSortedValA, const void * csrSortedRowPtrA, const void * csrSortedColIndA, const void * B, int ldb, csrsm2Info_t info, cusparseSolvePolicy_t policy, size_t * pBufferSize)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsm_analysis(cusparseHandle_t handle, int algo, cusparseOperation_t transA, cusparseOperation_t transB, int m, int nrhs, int nnz, const void * alpha, const cusparseMatDescr_t descrA, const void * csrSortedValA, const void * csrSortedRowPtrA, const void * csrSortedColIndA, const void * B, int ldb, csrsm2Info_t info, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsm_solve(cusparseHandle_t handle, int algo, cusparseOperation_t transA, cusparseOperation_t transB, int m, int nrhs, int nnz, const void * alpha, const cusparseMatDescr_t descrA, const void * csrSortedValA, const void * csrSortedRowPtrA, const void * csrSortedColIndA, const void * B, int ldb, csrsm2Info_t info, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (U*)alpha, descrA, (U*)csrSortedValA, (I*)csrSortedRowPtrA, (I*)csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
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
        cusparseSpMatDescr_t d_new;
        cusparseMatDescr_t d_leg;
        void * rowptrs = nullptr;
        void * colidxs = nullptr;
        void * vals = nullptr;
        int64_t nrows;
        int64_t ncols;
        int64_t nnz;
    };

    struct _descr_matrix_dense
    {
        cusparseDnMatDescr_t d_new;
        cusparseDnMatDescr_t d_new_complementary;
        void * vals = nullptr;
        int64_t nrows;
        int64_t ncols;
        int64_t ld;
        char order;
        descr_matrix_dense get_complementary()
        {
            descr_matrix_dense ret = std::make_shared<_descr_matrix_dense>();
            ret->d_new = d_new_complementary;
            ret->d_new_complementary = d_new;
            ret->vals = vals;
            ret->nrows = ncols;
            ret->ncols = nrows;
            ret->ld = ld;
            ret->order = mgm::change_order(order);
            return ret;
        }
    };

    struct _descr_vector_dense
    {
        cusparseDnVecDescr_t d_new;
        void * vals;
    };

    struct _descr_sparse_trsv
    {
        csrsv2Info_t i;
    };

    struct _descr_sparse_trsm
    {
        csrsm2Info_t i;
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
        CHECK(cusparseCreateCsr(&descr->d_new, nrows, ncols, nnz, dummyptr, dummyptr, dummyptr, _sparse_index_type<I>(), _sparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, _sparse_data_type<T>()));
        auto upper = CUSPARSE_FILL_MODE_UPPER;
        auto lower = CUSPARSE_FILL_MODE_LOWER;
        auto nonunit = CUSPARSE_DIAG_TYPE_NON_UNIT;
        if(symmetry == 'L') CHECK(cusparseSpMatSetAttribute(descr->d_new, CUSPARSE_SPMAT_FILL_MODE, &lower, sizeof(lower)));
        if(symmetry == 'U') CHECK(cusparseSpMatSetAttribute(descr->d_new, CUSPARSE_SPMAT_FILL_MODE, &upper, sizeof(upper)));
        if(symmetry != 'N') CHECK(cusparseSpMatSetAttribute(descr->d_new, CUSPARSE_SPMAT_DIAG_TYPE, &nonunit, sizeof(nonunit)));

        CHECK(cusparseCreateMatDescr(&descr->d_leg));
        CHECK(cusparseSetMatDiagType(descr->d_leg, CUSPARSE_DIAG_TYPE_NON_UNIT));
        CHECK(cusparseSetMatIndexBase(descr->d_leg, CUSPARSE_INDEX_BASE_ZERO));
        CHECK(cusparseSetMatType(descr->d_leg, CUSPARSE_MATRIX_TYPE_GENERAL));
        if(symmetry == 'U') CHECK(cusparseSetMatFillMode(descr->d_leg, CUSPARSE_FILL_MODE_UPPER));
        if(symmetry == 'L') CHECK(cusparseSetMatFillMode(descr->d_leg, CUSPARSE_FILL_MODE_LOWER));

        descr->nrows = nrows;
        descr->ncols = ncols;
        descr->nnz = nnz;
    }

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(cusparseCsrSetPointers(descr->d_new, matrix.rows, matrix.cols, matrix.vals));
        descr->rowptrs = matrix.rows;
        descr->colidxs = matrix.cols;
        descr->vals = matrix.vals;
    }

    void descr_matrix_csr_destroy(descr_matrix_csr & descr)
    {
        CHECK(cusparseDestroySpMat(descr->d_new));
        CHECK(cusparseDestroyMatDescr(descr->d_leg));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_matrix_dense_create(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order)
    {
        descr = std::make_shared<_descr_matrix_dense>();
        void * dummyptr = reinterpret_cast<void*>(1);
        if(order == 'R') CHECK(cusparseCreateDnMat(&descr->d_new, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_ROW));
        if(order == 'C') CHECK(cusparseCreateDnMat(&descr->d_new, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_COL));
        if(order == 'R') CHECK(cusparseCreateDnMat(&descr->d_new_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_COL));
        if(order == 'C') CHECK(cusparseCreateDnMat(&descr->d_new_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), CUSPARSE_ORDER_ROW));

        descr->nrows = nrows;
        descr->ncols = ncols;
        descr->ld = ld;
        descr->order = order;
    }

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(cusparseDnMatSetValues(descr->d_new, matrix.vals));
        CHECK(cusparseDnMatSetValues(descr->d_new_complementary, matrix.vals));
        descr->vals = matrix.vals;
    }

    void descr_matrix_dense_destroy(descr_matrix_dense & descr)
    {
        CHECK(cusparseDestroyDnMat(descr->d_new));
        CHECK(cusparseDestroyDnMat(descr->d_new_complementary));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_vector_dense_create(descr_vector_dense & descr, I size)
    {
        descr = std::make_shared<_descr_vector_dense>();
        void * dummyptr = reinterpret_cast<void*>(1);
        CHECK(cusparseCreateDnVec(&descr->d_new, size, dummyptr, _sparse_data_type<T>()));
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<T,I,A> & vector)
    {
        static_assert(A::is_data_device_accessible, "vector data must be device accessible");
        CHECK(cusparseDnVecSetValues(descr->d_new, vector.vals));
        descr->vals = vector.vals;
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx)
    {
        static_assert(A::is_data_device_accessible, "vector data must be device accessible");
        CHECK(cusparseDnVecSetValues(descr->d_new, matrix.vals + colidx * matrix.get_ld()));
        descr->vals = matrix.vals + colidx * matrix.get_ld();
    }

    void descr_vector_dense_destroy(descr_vector_dense & descr)
    {
        CHECK(cusparseDestroyDnVec(descr->d_new));
        descr.reset();
    }

    void descr_sparse_trsv_create(descr_sparse_trsv & descr)
    {
        descr = std::make_shared<_descr_sparse_trsv>();
        CHECK(cusparseCreateCsrsv2Info(&descr->i));
    }

    void descr_sparse_trsv_destroy(descr_sparse_trsv & descr)
    {
        CHECK(cusparseDestroyCsrsv2Info(descr->i));
        descr.reset();
    }

    void descr_sparse_trsm_create(descr_sparse_trsm & descr)
    {
        descr = std::make_shared<_descr_sparse_trsm>();
        CHECK(cusparseCreateCsrsm2Info(&descr->i));
    }

    void descr_sparse_trsm_destroy(descr_sparse_trsm & descr)
    {
        CHECK(cusparseDestroyCsrsm2Info(descr->i));
        descr.reset();
    }

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage)
    {
        if(transpose == 'T')
        {
            descr_matrix_dense descr_dense_complementary = dense->get_complementary();
            sparse_to_dense<T,I>(h, 'N', sparse, descr_dense_complementary, buffersize, buffer, stage);
            return;
        }
        if(stage == 'B') CHECK(cusparseSparseToDense_bufferSize(h->h, sparse->d_new, dense->d_new, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, &buffersize));
        if(stage == 'C') CHECK(cusparseSparseToDense           (h->h, sparse->d_new, dense->d_new, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, buffer));
    }

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage)
    {
        cusparseOperation_t op = (transpose == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        int bfsz = buffersize;
        T one = 1.0;
        if(stage == 'B') CHECK((_my_sparse_trsv_buffersize<T,I>)(h->h, op, matrix->nrows, matrix->nnz,       matrix->d_leg, matrix->vals, matrix->rowptrs, matrix->colidxs, descr_trsv->i, &bfsz));
        if(stage == 'P') CHECK((_my_sparse_trsv_analysis<T,I>)  (h->h, op, matrix->nrows, matrix->nnz,       matrix->d_leg, matrix->vals, matrix->rowptrs, matrix->colidxs, descr_trsv->i, policy, buffer));
        // if(stage == 'U') ;
        if(stage == 'C') CHECK((_my_sparse_trsv_solve<T,I>)     (h->h, op, matrix->nrows, matrix->nnz, &one, matrix->d_leg, matrix->vals, matrix->rowptrs, matrix->colidxs, descr_trsv->i, rhs->vals, sol->vals, policy, buffer));
        buffersize = bfsz;
    }

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
    {
        // cudaleg has the transB on the rhs as well as on the solution
        if(rhs->order != sol->order) eslog::error("dense matrix order has to be the same");
        if(rhs->order == 'R')
        {
            descr_matrix_dense descr_rhs_complementary = rhs->get_complementary();
            descr_matrix_dense descr_sol_complementary = sol->get_complementary();
            char transpose_rhs_compl = mgm::change_operation(transpose_rhs);
            trsm<T,I>(h, transpose_mat, transpose_rhs_compl, matrix, descr_rhs_complementary, descr_sol_complementary, descr_trsm, buffersize, buffer, stage);
            return;
        }
        cusparseOperation_t op_mat = (transpose_mat == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        cusparseOperation_t op_rhs = (transpose_rhs == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        T one = 1.0;
        int algo = 1;
        cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        I nrhs = (transpose_rhs == 'T' ? sol->nrows : sol->ncols);
        if(stage == 'B') CHECK((_my_sparse_trsm_buffersize<T,I>)(h->h, algo, op_mat, op_rhs, matrix->nrows, nrhs, matrix->nnz, &one, matrix->d_leg, matrix->vals, matrix->rowptrs, matrix->colidxs, rhs->vals, rhs->ld, descr_trsm->i, policy, &buffersize));
        if(stage == 'P') CHECK((_my_sparse_trsm_analysis<T,I>)  (h->h, algo, op_mat, op_rhs, matrix->nrows, nrhs, matrix->nnz, &one, matrix->d_leg, matrix->vals, matrix->rowptrs, matrix->colidxs, rhs->vals, rhs->ld, descr_trsm->i, policy, buffer));
        // if(stage == 'U') ;
        if(stage == 'C') CHECK((_my_sparse_trsm_solve<T,I>)     (h->h, algo, op_mat, op_rhs, matrix->nrows, nrhs, matrix->nnz, &one, matrix->d_leg, matrix->vals, matrix->rowptrs, matrix->colidxs, rhs->vals, rhs->ld, descr_trsm->i, policy, buffer));
        if(stage == 'C') CHECK(cudaMemcpy2DAsync(sol->vals, sol->ld * sizeof(T), rhs->vals, rhs->ld * sizeof(T), sol->nrows * sizeof(T), sol->ncols, cudaMemcpyDeviceToDevice, h->get_stream()));
    }

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage)
    {
        cusparseOperation_t op = (transpose == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        T one = 1.0;
        T zero = 0.0;
        if(stage == 'B') CHECK(cusparseSpMV_bufferSize(h->h, op, &one, A->d_new, x->d_new, &zero, y->d_new, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, &buffersize));
        if(stage == 'C') CHECK(cusparseSpMV           (h->h, op, &one, A->d_new, x->d_new, &zero, y->d_new, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, buffer));
    }

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        cusparseOperation_t op_A = (transpose_A == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        cusparseOperation_t op_B = (transpose_B == 'T' ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpMM_bufferSize(h->h, op_A, op_B, &one, A->d_new, B->d_new, &zero, C->d_new, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, &buffersize));
        if(stage == 'P') CHECK(cusparseSpMM_preprocess(h->h, op_A, op_B, &one, A->d_new, B->d_new, &zero, C->d_new, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
        if(stage == 'C') CHECK(cusparseSpMM           (h->h, op_A, op_B, &one, A->d_new, B->d_new, &zero, C->d_new, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
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

#ifdef MY_COMPILER_GCC
#pragma GCC diagnostic pop
#endif
#ifdef MY_COMPILER_CLANG
#pragma clang diagnostic pop
#endif

#undef MY_COMPILER_GCC
#undef MY_COMPILER_CLANG

#endif
#endif
