
#ifndef SRC_WRAPPERS_CUDA_W_CUDA_GPU_DNBLAS_HPP_
#define SRC_WRAPPERS_CUDA_W_CUDA_GPU_DNBLAS_HPP_

#ifdef HAVE_CUDA

#include "gpu/gpu_dnblas.h"

#include <cublas_v2.h>
#include "w.cuda.common.h"



namespace espreso {
namespace gpu {
namespace dnblas {

    namespace
    {
        inline static void _check(cublasStatus_t status, const char *file, int line)
        {
            if (status != CUBLAS_STATUS_SUCCESS)
            {
                char str[1000];
                snprintf(str, sizeof(str), "CUBLAS Error %d %s: %s. In file '%s' on line %d\n", status, cublasGetStatusName(status), cublasGetStatusString(status), file, line);
                eslog::error(str);
            }
        }

        template<typename T>
        static cublasStatus_t _my_blas_xtrsv(cublasHandle_t handle, cublasFillMode_t uplo, cublasOperation_t transA, cublasDiagType_t diag, int m, const T *A, int lda, T *x, int incx)
        {
            if constexpr(std::is_same_v<T,float>)  return cublasStrsv(handle, uplo, transA, diag, m, A, lda, x, incx);
            if constexpr(std::is_same_v<T,double>) return cublasDtrsv(handle, uplo, transA, diag, m, A, lda, x, incx);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xtrsm(cublasHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t transA, cublasDiagType_t diag, int m, int n, const T *alpha, const T *A, int lda, T *B, int ldb)
        {
            if constexpr(std::is_same_v<T,float>)  return cublasStrsm(handle, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
            if constexpr(std::is_same_v<T,double>) return cublasDtrsm(handle, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xsyrk(cublasHandle_t handle, cublasFillMode_t uplo, cublasOperation_t transA, int n, int k, const T *alpha, const T *A, int lda, const T *beta, T *C, int ldc)
        {
            if constexpr(std::is_same_v<T,float>)  return cublasSsyrk(handle, uplo, transA, n, k, alpha, A, lda, beta, C, ldc);
            if constexpr(std::is_same_v<T,double>) return cublasDsyrk(handle, uplo, transA, n, k, alpha, A, lda, beta, C, ldc);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xsymv(cublasHandle_t handle, cublasFillMode_t uplo, int n, const T *alpha, const T *A, int lda, const T *x, int incx, const T *beta, T *y, int incy)
        {
            if constexpr(std::is_same_v<T,float>)  return cublasSsymv(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
            if constexpr(std::is_same_v<T,double>) return cublasDsymv(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
        }
    }

    struct handle
    {
        cublasHandle_t h;
    };

    static void handle_create(handle & h, mgm::queue & q)
    {
        CHECK(cublasCreate(&h.h));
        CHECK(cublasSetStream(h.h, q.stream));
    }

    static void handle_destroy(handle & h)
    {
        CHECK(cublasDestroy(h.h));
    }

    template<typename F>
    void buffer_collect_size(handle & /*q*/, size_t & buffersize, const F & /*f*/)
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
    void trsv(handle & h, char mat_symmetry, char transpose, I n, I ld, T * matrix, T * rhs_sol)
    {
        cublasFillMode_t fill = (mat_symmetry == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        cublasOperation_t op = (transpose == 'T' ? CUBLAS_OP_T : CUBLAS_OP_N);
        CHECK(_my_blas_xtrsv(h.h, fill, op, CUBLAS_DIAG_NON_UNIT, n, matrix, ld, rhs_sol, 1));
    }

    template<typename T, typename I>
    void trsm(handle & h, char side, char mat_symmetry, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X)
    {
        cublasSideMode_t s = (side == 'L' ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT);
        cublasFillMode_t fill = (mat_symmetry == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        cublasOperation_t op = (transpose == 'T' ? CUBLAS_OP_T : CUBLAS_OP_N);
        T one = 1.0;
        CHECK(_my_blas_xtrsm(h.h, s, fill, op, CUBLAS_DIAG_NON_UNIT, nrows_X, ncols_X, &one, A, ld_A, rhs_sol, ld_X));
    }

    template<typename T, typename I>
    void syrk(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C)
    {
        cublasFillMode_t fill = (out_fill == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        cublasOperation_t op = (transpose == 'T' ? CUBLAS_OP_T : CUBLAS_OP_N);
        T zero = 0.0;
        T one = 1.0;
        CHECK(_my_blas_xsyrk(h.h, fill, op, n, k, &one, A, ld_A, &zero, C, ld_C));
    }

    template<typename T, typename I>
    void symv(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out)
    {
        cublasFillMode_t f = (fill == 'U' ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
        T zero = 0.0;
        T one = 1.0;
        CHECK(_my_blas_xsymv(h.h, f, n, &one, A, ld_A, vec_in, 1, &zero, vec_out, 1));
    }

}
}
}

#endif
#endif /* SRC_WRAPPERS_CUDA_W_CUDA_GPU_DNBLAS_HPP_ */
