
#ifdef HAVE_CUDA

#include "gpu/gpu_dnblas.h"
#include "w.cuda.gpu_management.h"
#include "basis/utilities/utils.h"

#include <cublas_v2.h>
#include <complex>



inline void _check(cublasStatus_t status, const char *file, int line)
{
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        char str[1000];
        snprintf(str, sizeof(str), "CUBLAS Error %d %s: %s. In file '%s' on line %d\n", status, cublasGetStatusName(status), cublasGetStatusString(status), file, line);
        eslog::error(str);
    }
}



namespace espreso {
namespace gpu {
namespace dnblas {

    namespace
    {
        template<typename T> struct cpp_to_cuda_type { using type = T; };
        template<> struct cpp_to_cuda_type<std::complex<float>> { using type = cuComplex; };
        template<> struct cpp_to_cuda_type<std::complex<double>> { using type = cuDoubleComplex; };
        template<typename T> using cpp_to_cuda_type_t = typename cpp_to_cuda_type<T>::type;

        template<typename T>
        static cublasStatus_t _my_blas_xtrsv(cublasHandle_t handle, cublasFillMode_t uplo, cublasOperation_t transA, cublasDiagType_t diag, int m, const T *A, int lda, T *x, int incx)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return cublasStrsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
            if constexpr(std::is_same_v<T,double>)               return cublasDtrsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return cublasCtrsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<double>>) return cublasZtrsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xtrsm(cublasHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t transA, cublasDiagType_t diag, int m, int n, const T *alpha, const T *A, int lda, T *B, int ldb)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return cublasStrsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
            if constexpr(std::is_same_v<T,double>)               return cublasDtrsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return cublasCtrsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
            if constexpr(std::is_same_v<T,std::complex<double>>) return cublasZtrsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xsyrk(cublasHandle_t handle, cublasFillMode_t uplo, cublasOperation_t transA, int n, int k, const T *alpha, const T *A, int lda, const T *beta, T *C, int ldc)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,float>)  return cublasSsyrk(handle, uplo, transA, n, k, (U*)alpha, (U*)A, lda, (U*)beta, (U*)C, ldc);
            if constexpr(std::is_same_v<T,double>) return cublasDsyrk(handle, uplo, transA, n, k, (U*)alpha, (U*)A, lda, (U*)beta, (U*)C, ldc);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xherk(cublasHandle_t handle, cublasFillMode_t uplo, cublasOperation_t transA, int n, int k, const utils::remove_complex_t<T> *alpha, const T *A, int lda, const utils::remove_complex_t<T> *beta, T *C, int ldc)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,std::complex<float>>)  return cublasCherk(handle, uplo, transA, n, k, alpha, (U*)A, lda, beta, (U*)C, ldc);
            if constexpr(std::is_same_v<T,std::complex<double>>) return cublasZherk(handle, uplo, transA, n, k, alpha, (U*)A, lda, beta, (U*)C, ldc);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xsymv(cublasHandle_t handle, cublasFillMode_t uplo, int n, const T *alpha, const T *A, int lda, const T *x, int incx, const T *beta, T *y, int incy)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,float>)  return cublasSsymv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,double>) return cublasDsymv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xhemv(cublasHandle_t handle, cublasFillMode_t uplo, int n, const T *alpha, const T *A, int lda, const T *x, int incx, const T *beta, T *y, int incy)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,std::complex<float>>)  return cublasChemv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,std::complex<double>>) return cublasZhemv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
        }
    }

    struct _handle
    {
        cublasHandle_t h;
    };

    static void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_unique<_handle>();
        CHECK(cublasCreate(&h->h));
        CHECK(cublasSetStream(h->h, q->stream));
    }

    static void handle_destroy(handle & h)
    {
        CHECK(cublasDestroy(h->h));
        h.reset();
    }

    static void buffer_collect_size(handle & /*q*/, size_t & buffersize, const std::function<void(void)> & /*f*/)
    {
        // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
        // no manual workspace needed if I use just a single stream with this handle
        buffersize = 0;
    }

    static void buffer_set(handle & /*q*/, void * /*ptr*/, size_t /*size*/)
    {
    }

    static void buffer_unset(handle & /*q*/)
    {
    }

    template<typename T, typename I>
    static void trsv(handle & h, char mat_symmetry, char transpose, I n, I ld, T * matrix, T * rhs_sol)
    {
        cublasFillMode_t fill = (mat_symmetry == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        cublasOperation_t op = (transpose == 'T' ? CUBLAS_OP_T : CUBLAS_OP_N);
        CHECK(_my_blas_xtrsv<T>(h->h, fill, op, CUBLAS_DIAG_NON_UNIT, n, matrix, ld, rhs_sol, 1));
    }

    template<typename T, typename I>
    static void trsm(handle & h, char side, char mat_symmetry, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X)
    {
        cublasSideMode_t s = (side == 'L' ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT);
        cublasFillMode_t fill = (mat_symmetry == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        cublasOperation_t op = (transpose == 'T' ? CUBLAS_OP_T : CUBLAS_OP_N);
        T one = 1.0;
        CHECK(_my_blas_xtrsm<T>(h->h, s, fill, op, CUBLAS_DIAG_NON_UNIT, nrows_X, ncols_X, &one, A, ld_A, rhs_sol, ld_X));
    }

    template<typename T, typename I>
    static void herk(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C)
    {
        cublasFillMode_t fill = (out_fill == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        cublasOperation_t op = (transpose == 'T' ? CUBLAS_OP_T : CUBLAS_OP_N);
        utils::remove_complex_t<T> zero = 0.0;
        utils::remove_complex_t<T> one = 1.0;
        if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsyrk<T>(h->h, fill, op, n, k, &one, A, ld_A, &zero, C, ld_C));
        if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xherk<T>(h->h, fill, op, n, k, &one, A, ld_A, &zero, C, ld_C));
    }

    template<typename T, typename I>
    static void hemv(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out)
    {
        cublasFillMode_t f = (fill == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        T zero = 0.0;
        T one = 1.0;
        if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsymv<T>(h->h, f, n, &one, A, ld_A, vec_in, 1, &zero, vec_out, 1));
        if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xhemv<T>(h->h, f, n, &one, A, ld_A, vec_in, 1, &zero, vec_out, 1));
    }



    #define INSTANTIATE(T,I) \
    template void trsv<T,I>(handle & h, char mat_symmetry, char transpose, I n, I ld, T * matrix, T * rhs_sol); \
    template void trsm<T,I>(handle & h, char side, char mat_symmetry, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X); \
    template void herk<T,I>(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C); \
    template void hemv<T,I>(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out);
        // INSTANTIATE(float,                int32_t)
        INSTANTIATE(double,               int32_t)
        // INSTANTIATE(std::complex<float >, int32_t)
        // INSTANTIATE(std::complex<double>, int32_t)
        // INSTANTIATE(float,                int64_t)
        // INSTANTIATE(double,               int64_t)
        // INSTANTIATE(std::complex<float >, int64_t)
        // INSTANTIATE(std::complex<double>, int64_t)
    #undef INSTANTIATE

}
}
}

#endif
