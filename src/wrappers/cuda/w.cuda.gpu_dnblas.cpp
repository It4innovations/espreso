
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
        espreso::eslog::error("CUBLAS Error %d %s: %s. In file '%s' on line %d\n", status, cublasGetStatusName(status), cublasGetStatusString(status), file, line);
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

        static cublasFillMode_t _char_to_fill(char c)
        {
            switch(c)
            {
                case 'U': return CUBLAS_FILL_MODE_UPPER;
                case 'L': return CUBLAS_FILL_MODE_LOWER;
                default: eslog::error("invalid fill '%c'\n", c);
            }
        }

        static cublasOperation_t _char_to_operation(char c)
        {
            switch(c)
            {
                case 'N': return CUBLAS_OP_N;
                case 'T': return CUBLAS_OP_T;
                case 'H': return CUBLAS_OP_C;
                default: eslog::error("invalid operation '%c'\n", c);
            }
        }

        static cublasSideMode_t _char_to_side(char c)
        {
            switch(c)
            {
                case 'L': return CUBLAS_SIDE_LEFT;
                case 'R': return CUBLAS_SIDE_RIGHT;
                default: eslog::error("invalid side '%c'\n", c);
            }
        }
    }

    struct _handle
    {
        cublasHandle_t h;
    };

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        CHECK(cublasCreate(&h->h));
        CHECK(cublasSetStream(h->h, q->stream));
    }

    void handle_destroy(handle & h)
    {
        CHECK(cublasDestroy(h->h));
        h.reset();
    }

    void buffer_collect_size(handle & /*q*/, size_t & buffersize, const std::function<void(void)> & /*f*/)
    {
        // https://docs.nvidia.com/cuda/cublas/index.html#cublassetworkspace
        // no manual workspace needed if I use just a single stream with this handle
        buffersize = 0;
    }

    void buffer_set(handle & /*q*/, void * /*ptr*/, size_t /*size*/)
    {
    }

    void buffer_unset(handle & /*q*/)
    {
    }

    template<typename T, typename I>
    void trsv(handle & h, char fill, char transpose, I n, I ld, T * matrix, T * rhs_sol)
    {
        CHECK(_my_blas_xtrsv<T>(h->h, _char_to_fill(fill), _char_to_operation(transpose), CUBLAS_DIAG_NON_UNIT, n, matrix, ld, rhs_sol, 1));
    }

    template<typename T, typename I>
    void trsm(handle & h, char side, char fill, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X)
    {
        T one = 1.0;
        CHECK(_my_blas_xtrsm<T>(h->h, _char_to_side(side), _char_to_fill(fill), _char_to_operation(transpose), CUBLAS_DIAG_NON_UNIT, nrows_X, ncols_X, &one, A, ld_A, rhs_sol, ld_X));
    }

    template<typename T, typename I>
    void herk(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C)
    {
        if(utils::is_real<T>() && transpose == 'H') transpose = 'T';
        utils::remove_complex_t<T> zero = 0.0;
        utils::remove_complex_t<T> one = 1.0;
        if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsyrk<T>(h->h, _char_to_fill(out_fill), _char_to_operation(transpose), n, k, &one, A, ld_A, &zero, C, ld_C));
        if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xherk<T>(h->h, _char_to_fill(out_fill), _char_to_operation(transpose), n, k, &one, A, ld_A, &zero, C, ld_C));
    }

    template<typename T, typename I>
    void hemv(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out)
    {
        T zero = 0.0;
        T one = 1.0;
        if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsymv<T>(h->h, _char_to_fill(fill), n, &one, A, ld_A, vec_in, 1, &zero, vec_out, 1));
        if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xhemv<T>(h->h, _char_to_fill(fill), n, &one, A, ld_A, vec_in, 1, &zero, vec_out, 1));
    }



    #define INSTANTIATE_T_I(T,I) \
    template void trsv<T,I>(handle & h, char fill, char transpose, I n, I ld, T * matrix, T * rhs_sol); \
    template void trsm<T,I>(handle & h, char side, char fill, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X); \
    template void herk<T,I>(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C); \
    template void hemv<T,I>(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out);

        #define INSTANTIATE_T(T) \
        INSTANTIATE_T_I(T, int32_t) \
        /* INSTANTIATE_T_I(T, int64_t) */

            // INSTANTIATE_T(float)
            INSTANTIATE_T(double)
            // INSTANTIATE_T(std::complex<float>)
            // INSTANTIATE_T(std::complex<double>)

        #undef INSTANTIATE_T
    #undef INSTANTIATE_T_I

}
}
}

#endif
