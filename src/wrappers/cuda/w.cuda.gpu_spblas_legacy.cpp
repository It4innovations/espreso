
#ifdef HAVE_CUDA
#ifdef USE_CUSPARSE_LEGACY

#include "gpu/gpu_spblas.h"
#include "w.cuda.gpu_management.h"
#include "basis/utilities/utils.h"

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
        espreso::eslog::error("CUSPARSE Error %d %s: %s. In file '%s' on line %d\n", status, cusparseGetErrorName(status), cusparseGetErrorString(status), file, line);
    }
}



#include "w.cuda.gpu_spblas_common.h"



namespace espreso {
namespace gpu {
namespace spblas {

    spblas_wrapper_impl get_implementation()
    {
        return spblas_wrapper_impl::CUSPARSE_LEGACY;
    }

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
        static cusparseStatus_t _my_sparse_trsv_buffersize(cusparseHandle_t handle, cusparseOperation_t transA, I m, I nnz, const cusparseMatDescr_t descrA, T * csrValA, const I * csrRowPtrA, const I * csrColIndA, csrsv2Info_t info, int * pBufferSizeInBytes)
        {
            static_assert(std::is_same_v<I,int32_t>, "only 32-bit integers are supported in legacy cusparse");
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, pBufferSizeInBytes);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, pBufferSizeInBytes);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, pBufferSizeInBytes);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsv2_bufferSize(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, pBufferSizeInBytes);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsv_analysis(cusparseHandle_t handle, cusparseOperation_t transA, I m, I nnz, const cusparseMatDescr_t descrA, T * csrValA, const I * csrRowPtrA, const I * csrColIndA, csrsv2Info_t info, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            static_assert(std::is_same_v<I,int32_t>, "only 32-bit integers are supported in legacy cusparse");
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsv2_analysis(handle, transA, m, nnz, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, policy, pBuffer);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsv_solve(cusparseHandle_t handle, cusparseOperation_t transA, I m, I nnz, const T * alpha, const cusparseMatDescr_t descrA, T * csrValA, const I * csrRowPtrA, const I * csrColIndA, csrsv2Info_t info, const T * x, T * y, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            static_assert(std::is_same_v<I,int32_t>, "only 32-bit integers are supported in legacy cusparse");
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsv2_solve(handle, transA, m, nnz, (const U*)alpha, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, (const U*)x, (U*)y, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsv2_solve(handle, transA, m, nnz, (const U*)alpha, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, (const U*)x, (U*)y, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsv2_solve(handle, transA, m, nnz, (const U*)alpha, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, (const U*)x, (U*)y, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsv2_solve(handle, transA, m, nnz, (const U*)alpha, descrA, (U*)csrValA, csrRowPtrA, csrColIndA, info, (const U*)x, (U*)y, policy, pBuffer);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsm_buffersize(cusparseHandle_t handle, int algo, cusparseOperation_t transA, cusparseOperation_t transB, int m, int nrhs, int nnz, const T * alpha, const cusparseMatDescr_t descrA, const T * csrSortedValA, const I * csrSortedRowPtrA, const I * csrSortedColIndA, T * B, int ldb, csrsm2Info_t info, cusparseSolvePolicy_t policy, size_t * pBufferSize)
        {
            static_assert(std::is_same_v<I,int32_t>, "only 32-bit integers are supported in legacy cusparse");
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsm2_bufferSizeExt(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBufferSize);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsm_analysis(cusparseHandle_t handle, int algo, cusparseOperation_t transA, cusparseOperation_t transB, int m, int nrhs, int nnz, const T * alpha, const cusparseMatDescr_t descrA, const T * csrSortedValA, const I * csrSortedRowPtrA, const I * csrSortedColIndA, T * B, int ldb, csrsm2Info_t info, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            static_assert(std::is_same_v<I,int32_t>, "only 32-bit integers are supported in legacy cusparse");
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsm2_analysis(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
        }

        template<typename T, typename I>
        static cusparseStatus_t _my_sparse_trsm_solve(cusparseHandle_t handle, int algo, cusparseOperation_t transA, cusparseOperation_t transB, int m, int nrhs, int nnz, const T * alpha, const cusparseMatDescr_t descrA, const T * csrSortedValA, const I * csrSortedRowPtrA, const I * csrSortedColIndA, T * B, int ldb, csrsm2Info_t info, cusparseSolvePolicy_t policy, void * pBuffer)
        {
            static_assert(std::is_same_v<I,int32_t>, "only 32-bit integers are supported in legacy cusparse");
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T, float>)                return cusparseScsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, double>)               return cusparseDcsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<float>>)  return cusparseCcsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
            if constexpr(std::is_same_v<T, std::complex<double>>) return cusparseZcsrsm2_solve(handle, algo, transA, transB, m, nrhs, nnz, (const U*)alpha, descrA, (const U*)csrSortedValA, csrSortedRowPtrA, csrSortedColIndA, (U*)B, ldb, info, policy, pBuffer);
        }

        template<typename T>
        static cusparseOperation_t _char_to_operation(char c)
        {
            if(utils::is_real<T>())
            {
                switch(c)
                {
                    case 'N': return CUSPARSE_OPERATION_NON_TRANSPOSE;
                    case 'T': return CUSPARSE_OPERATION_TRANSPOSE;
                    case 'H': return CUSPARSE_OPERATION_TRANSPOSE;
                    default: eslog::error("invalid operation '%c'\n", c);
                }
            }
            if(utils::is_complex<T>())
            {
                switch(c)
                {
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
        _descr_matrix_dense get_complementary()
        {
            _descr_matrix_dense ret;
            ret.d_new = d_new_complementary;
            ret.d_new_complementary = d_new;
            ret.vals = vals;
            ret.nrows = ncols;
            ret.ncols = nrows;
            ret.ld = ld;
            ret.order = mgm::order_change(order);
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
        if(h.get() == nullptr) return;

        CHECK(cusparseDestroy(h->h));
        h.reset();
    }

    template<typename T, typename I>
    void descr_matrix_csr_create(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill)
    {
        descr = std::make_shared<_descr_matrix_csr>();
        void * dummyptr = reinterpret_cast<void*>(1);
        CHECK(cusparseCreateCsr(&descr->d_new, nrows, ncols, nnz, dummyptr, dummyptr, dummyptr, _sparse_index_type<I>(), _sparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, _sparse_data_type<T>()));
        auto upper = CUSPARSE_FILL_MODE_UPPER;
        auto lower = CUSPARSE_FILL_MODE_LOWER;
        auto nonunit = CUSPARSE_DIAG_TYPE_NON_UNIT;
        if(fill == 'L') CHECK(cusparseSpMatSetAttribute(descr->d_new, CUSPARSE_SPMAT_FILL_MODE, &lower, sizeof(lower)));
        if(fill == 'U') CHECK(cusparseSpMatSetAttribute(descr->d_new, CUSPARSE_SPMAT_FILL_MODE, &upper, sizeof(upper)));
        if(fill != 'N') CHECK(cusparseSpMatSetAttribute(descr->d_new, CUSPARSE_SPMAT_DIAG_TYPE, &nonunit, sizeof(nonunit)));

        CHECK(cusparseCreateMatDescr(&descr->d_leg));
        CHECK(cusparseSetMatDiagType(descr->d_leg, CUSPARSE_DIAG_TYPE_NON_UNIT));
        CHECK(cusparseSetMatIndexBase(descr->d_leg, CUSPARSE_INDEX_BASE_ZERO));
        CHECK(cusparseSetMatType(descr->d_leg, CUSPARSE_MATRIX_TYPE_GENERAL));
        if(fill == 'U') CHECK(cusparseSetMatFillMode(descr->d_leg, CUSPARSE_FILL_MODE_UPPER));
        if(fill == 'L') CHECK(cusparseSetMatFillMode(descr->d_leg, CUSPARSE_FILL_MODE_LOWER));

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
        if(descr.get() == nullptr) return;

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
        if(descr.get() == nullptr) return;

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
        if(descr.get() == nullptr) return;

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
        if(descr.get() == nullptr) return;

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
        if(descr.get() == nullptr) return;

        CHECK(cusparseDestroyCsrsm2Info(descr->i));
        descr.reset();
    }

    template<typename T, typename I>
    void transpose(handle & h, descr_matrix_csr & output, descr_matrix_csr & input, bool conjugate, size_t & buffersize, void * buffer, char stage)
    {
        cudaStream_t stream = h->get_stream();
        if(stage == 'B') my_csr_transpose_buffersize<I>(stream, input->nrows, input->ncols, input->nnz, buffersize);
        if(stage == 'P') my_csr_transpose_preprocess<I>(stream, input->nrows, input->ncols, input->nnz, (I*)input->rowptrs, (I*)input->colidxs, (I*)output->rowptrs, (I*)output->colidxs, buffersize, buffer);
        if(stage == 'C') my_csr_transpose_compute<T,I>(stream, input->nnz, (T*)input->vals, (T*)output->vals, conjugate, buffer);
    }

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage)
    {
        if(transpose == 'N')
        {
            if(stage == 'B') CHECK(cusparseSparseToDense_bufferSize(h->h, sparse->d_new, dense->d_new, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, &buffersize));
            if(stage == 'C') CHECK(cusparseSparseToDense           (h->h, sparse->d_new, dense->d_new, CUSPARSE_SPARSETODENSE_ALG_DEFAULT, buffer));
        }
        else if(transpose == 'T')
        {
            descr_matrix_dense descr_dense_complementary = std::make_shared<_descr_matrix_dense>(dense->get_complementary());
            sparse_to_dense<T,I>(h, 'N', sparse, descr_dense_complementary, buffersize, buffer, stage);
        }
        else eslog::error("transpose '%c' not supported\n", transpose);
    }

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage)
    {
        cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        int bfsz = buffersize;
        T one = 1.0;
        if(stage == 'B') CHECK((_my_sparse_trsv_buffersize<T,I>)(h->h, _char_to_operation<T>(transpose), matrix->nrows, matrix->nnz,       matrix->d_leg, (T*)matrix->vals, (I*)matrix->rowptrs, (I*)matrix->colidxs, descr_trsv->i, &bfsz));
        if(stage == 'P') CHECK((_my_sparse_trsv_analysis<T,I>)  (h->h, _char_to_operation<T>(transpose), matrix->nrows, matrix->nnz,       matrix->d_leg, (T*)matrix->vals, (I*)matrix->rowptrs, (I*)matrix->colidxs, descr_trsv->i, policy, buffer));
        // if(stage == 'U') ;
        if(stage == 'C') CHECK((_my_sparse_trsv_solve<T,I>)     (h->h, _char_to_operation<T>(transpose), matrix->nrows, matrix->nnz, &one, matrix->d_leg, (T*)matrix->vals, (I*)matrix->rowptrs, (I*)matrix->colidxs, descr_trsv->i, (T*)rhs->vals, (T*)sol->vals, policy, buffer));
        buffersize = bfsz;
    }

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
    {
        // legacy cusparse assumes colmajor for both dense matrices
        // legacy cusparse api has the transB on the rhs as well as on the solution
        if(rhs->order == sol->order)
        {
            if(transpose_rhs == transpose_sol)
            {
                if(rhs->order == 'C')
                {
                    T one = 1.0;
                    int algo = 1;
                    cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
                    I nrhs = (transpose_rhs == 'T' ? sol->nrows : sol->ncols);
                    if(stage == 'B') CHECK((_my_sparse_trsm_buffersize<T,I>)(h->h, algo, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), matrix->nrows, nrhs, matrix->nnz, &one, matrix->d_leg, (T*)matrix->vals, (I*)matrix->rowptrs, (I*)matrix->colidxs, (T*)rhs->vals, rhs->ld, descr_trsm->i, policy, &buffersize));
                    if(stage == 'P') CHECK((_my_sparse_trsm_analysis<T,I>)  (h->h, algo, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), matrix->nrows, nrhs, matrix->nnz, &one, matrix->d_leg, (T*)matrix->vals, (I*)matrix->rowptrs, (I*)matrix->colidxs, (T*)rhs->vals, rhs->ld, descr_trsm->i, policy, buffer));
                    // if(stage == 'U') ; // I can just modify the values behind the pointer, and cusparse will notice
                    if(stage == 'C') CHECK((_my_sparse_trsm_solve<T,I>)     (h->h, algo, _char_to_operation<T>(transpose_mat), _char_to_operation<T>(transpose_rhs), matrix->nrows, nrhs, matrix->nnz, &one, matrix->d_leg, (T*)matrix->vals, (I*)matrix->rowptrs, (I*)matrix->colidxs, (T*)rhs->vals, rhs->ld, descr_trsm->i, policy, buffer));
                    if(stage == 'C') CHECK(cudaMemcpy2DAsync(sol->vals, sol->ld * sizeof(T), rhs->vals, rhs->ld * sizeof(T), sol->nrows * sizeof(T), sol->ncols, cudaMemcpyDeviceToDevice, h->get_stream()));
                }
                else if(rhs->order == 'R')
                {
                    descr_matrix_dense descr_rhs_compl = std::make_shared<_descr_matrix_dense>(rhs->get_complementary());
                    descr_matrix_dense descr_sol_compl = std::make_shared<_descr_matrix_dense>(sol->get_complementary());
                    char transpose_compl = mgm::operation_combine(transpose_rhs, 'T');
                    trsm<T,I>(h, transpose_mat, transpose_compl, transpose_compl, matrix, descr_rhs_compl, descr_sol_compl, descr_trsm, buffersize, buffer, stage);
                }
                else
                {
                    eslog::error("wrong dense matrix order '%c'\n", rhs->order);
                }
            }
            else
            {
                eslog::error("unsupported combimation of matrix orders and transpositions '%c'\n", rhs->order);
            }
        }
        else
        {
            descr_matrix_dense descr_rhs_compl = std::make_shared<_descr_matrix_dense>(rhs->get_complementary());
            char transpose_rhs_compl = mgm::operation_combine(transpose_rhs, 'T');
            trsm<T,I>(h, transpose_mat, transpose_rhs_compl, transpose_sol, matrix, descr_rhs_compl, sol, descr_trsm, buffersize, buffer, stage);
        }
    }

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage)
    {
        T one = 1.0;
        T zero = 0.0;
        if(stage == 'B') CHECK(cusparseSpMV_bufferSize(h->h, _char_to_operation<T>(transpose), &one, A->d_new, x->d_new, &zero, y->d_new, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, &buffersize));
        if(stage == 'C') CHECK(cusparseSpMV           (h->h, _char_to_operation<T>(transpose), &one, A->d_new, x->d_new, &zero, y->d_new, _sparse_data_type<T>(), CUSPARSE_SPMV_ALG_DEFAULT, buffer));
    }

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'B') CHECK(cusparseSpMM_bufferSize(h->h, _char_to_operation<T>(transpose_A), _char_to_operation<T>(transpose_B), &one, A->d_new, B->d_new, &zero, C->d_new, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, &buffersize));
        if(stage == 'P') CHECK(cusparseSpMM_preprocess(h->h, _char_to_operation<T>(transpose_A), _char_to_operation<T>(transpose_B), &one, A->d_new, B->d_new, &zero, C->d_new, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
        if(stage == 'C') CHECK(cusparseSpMM           (h->h, _char_to_operation<T>(transpose_A), _char_to_operation<T>(transpose_B), &one, A->d_new, B->d_new, &zero, C->d_new, _sparse_data_type<T>(), CUSPARSE_SPMM_ALG_DEFAULT, buffer));
    }

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

#include "gpu/gpu_spblas.inst.hpp"

#endif
#endif
