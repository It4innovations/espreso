
#ifdef HAVE_CUDA
#ifdef ESPRESO_USE_WRAPPER_GPU_CUDA

#include "gpu/gpu_dnblas.h"
#include "w.cuda.gpu_management.h"
#include "basis/utilities/utils.h"

#include <cublas_v2.h>
#include <complex>



inline void _check(cublasStatus_t status, const char *file, int line)
{
    if (status != CUBLAS_STATUS_SUCCESS) {
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

        template<typename T>
        static cublasStatus_t _my_blas_xgemv(cublasHandle_t handle, cublasOperation_t trans, int m, int n, const T *alpha, const T *A, int lda, const T *x, int incx, const T *beta, T *y, int incy)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return cublasSgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,double>)               return cublasDgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return cublasCgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,std::complex<double>>) return cublasZgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
        }

        template<typename T>
        static cublasStatus_t _my_blas_xscal(cublasHandle_t handle, int n, const T *alpha, T *x, int incx)
        {
            using U = cpp_to_cuda_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return cublasSscal(handle, n, (U*)alpha, (U*)x, incx);
            if constexpr(std::is_same_v<T,double>)               return cublasDscal(handle, n, (U*)alpha, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return cublasCscal(handle, n, (U*)alpha, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<double>>) return cublasZscal(handle, n, (U*)alpha, (U*)x, incx);
        }

        static cublasFillMode_t _char_to_fill(char c)
        {
            switch(c) {
                case 'U': return CUBLAS_FILL_MODE_UPPER;
                case 'L': return CUBLAS_FILL_MODE_LOWER;
                default: eslog::error("invalid fill '%c'\n", c);
            }
        }

        static cublasOperation_t _char_to_operation(char c)
        {
            switch(c) {
                case 'N': return CUBLAS_OP_N;
                case 'T': return CUBLAS_OP_T;
                case 'H': return CUBLAS_OP_C;
                default: eslog::error("invalid operation '%c'\n", c);
            }
        }

        static cublasSideMode_t _char_to_side(char c)
        {
            switch(c) {
                case 'L': return CUBLAS_SIDE_LEFT;
                case 'R': return CUBLAS_SIDE_RIGHT;
                default: eslog::error("invalid side '%c'\n", c);
            }
        }
    }

    void init_library(mgm::queue & q) {}

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        CHECK(cublasCreate(&h->h));
        CHECK(cublasSetStream(h->h, q->stream));
    }

    void handle_destroy(handle & h)
    {
        if(h.get() == nullptr) return;

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
    void trsv(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x)
    {
        if(order_A == 'C') {
            if(op_A == 'C') { // conjugate_only operation not supported ...
                utils::remove_complex_t<T> neg_one = -1;
                if constexpr(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(x) + 1, 2);
                trsv<T,I>(h, n, A, ld_A, order_A, 'N', fill_A, x);
                if constexpr(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(x) + 1, 2);
            }
            else {
                CHECK(_my_blas_xtrsv(h->h, _char_to_fill(fill_A), _char_to_operation(op_A), CUBLAS_DIAG_NON_UNIT, n, A, ld_A, x, 1));
            }
        }
        else if(order_A == 'R') {
            char op_A_compl = mgm::operation_combine(op_A, 'T');
            char fill_A_compl = mgm::fill_change(fill_A);
            trsv<T,I>(h, n, A, ld_A, 'C', op_A_compl, fill_A_compl, x);
        }
        else {
            eslog::error("invalid order_A '%c'\n", order_A);
        }
    }

    template<typename T, typename I>
    void trsm(handle & h, char side, I n, I nrhs, T * A, I ld_A, char order_A, char op_A, char fill_A, T * X, I ld_X, char order_X, char op_X)
    {
        // eslog::info("TRSM side %c order_A %c op_A %c fill_A %c order_X %c op_X %c\n", side, order_A, op_A, fill_A, order_X, op_X);
        // solve op(A) op(X) = op(B) or op(X) op(A) = op(B), where input B is stored in X
        if(order_A == 'C') {
            if(order_X == 'C') {
                if(op_X == 'N') {
                    if(op_A == 'C') { // conjugate_only operation not supported ...
                        utils::remove_complex_t<T> neg_one = -1;
                        I nvals_X = (side == 'L' ? nrhs : n) * ld_X; // i know order_X == 'C'
                        if constexpr(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, nvals_X, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(X) + 1, 2);
                        trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, 'N', fill_A, X, ld_X, order_X, op_X);
                        if constexpr(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, nvals_X, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(X) + 1, 2);
                    }
                    else {
                        T one = 1.0;
                        I nrows_X = (side == 'L' ? n : nrhs);
                        I ncols_X = (side == 'L' ? nrhs : n);
                        CHECK(_my_blas_xtrsm<T>(h->h, _char_to_side(side), _char_to_fill(fill_A), _char_to_operation(op_A), CUBLAS_DIAG_NON_UNIT, nrows_X, ncols_X, &one, A, ld_A, X, ld_X));
                    }
                }
                else if(op_X == 'C') {
                    char op_A_comb = mgm::operation_combine(op_A, 'C');
                    trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, op_A_comb, fill_A, X, ld_X, order_X, 'N');
                }
                else if(op_X == 'T' || op_X == 'H') {
                    char side_compl = mgm::side_change(side);
                    char op_A_comb = mgm::operation_combine(op_A, op_X);
                    trsm<T,I>(h, side_compl, n, nrhs, A, ld_A, order_A, op_A_comb, fill_A, X, ld_X, order_X, 'N');
                }
                else {
                    eslog::error("invalid op_X %c\n", op_X);
                }
            }
            else if(order_X == 'R') {
                char op_X_compl = mgm::operation_combine(op_X, 'T');
                trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, op_A, fill_A, X, ld_X, 'C', op_X_compl);
            }
            else {
                eslog::error("invalid order_X '%c'\n", order_X);
            }
        }
        else if(order_A == 'R') {
            char op_A_compl = mgm::operation_combine(op_A, 'T');
            char fill_A_compl = mgm::fill_change(fill_A);
            trsm<T,I>(h, side, n, nrhs, A, ld_A, 'C', op_A_compl, fill_A_compl, X, ld_X, order_X, op_X);
        }
        else {
            eslog::error("invalid order_A '%c'\n", order_A);
        }
    }

    template<typename T, typename I>
    void herk(handle & h, I n, I k, T * A, I ld_A, char order_A, char op_A, T * C, I ld_C, char order_C, char fill_C)
    {
        // C = op(A) * op(A)^H
        if(order_C == 'C') {
            if(order_A == 'C') {
                if(utils::is_real<T>() && op_A == 'H') op_A = 'T';
                utils::remove_complex_t<T> zero = 0.0;
                utils::remove_complex_t<T> one = 1.0;
                if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsyrk<T>(h->h, _char_to_fill(fill_C), _char_to_operation(op_A), n, k, &one, A, ld_A, &zero, C, ld_C));
                if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xherk<T>(h->h, _char_to_fill(fill_C), _char_to_operation(op_A), n, k, &one, A, ld_A, &zero, C, ld_C));
            }
            else if(order_A == 'R') {
                char op_A_compl = mgm::operation_combine(op_A, 'T');
                herk<T,I>(h, n, k, A, ld_A, 'C', op_A_compl, C, ld_C, order_C, fill_C);
            }
            else {
                eslog::error("invalid order_A '%c'\n", order_A);
            }
        }
        else if(order_C == 'R') {
            char fill_C_compl = mgm::fill_change(fill_C);
            char op_A_compl = mgm::operation_combine(op_A, 'C');
            herk<T,I>(h, n, k, A, ld_A, order_A, op_A_compl, C, ld_C, 'C', fill_C_compl);
        }
        else {
            eslog::error("invalid order_C '%c'\n", order_C);
        }
    }

    template<typename T, typename I>
    void hemv(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x, T * y)
    {
        if(order_A == 'C') {
            if constexpr(utils::is_real<T>()) {
                T zero = 0.0;
                T one = 1.0;
                CHECK(_my_blas_xsymv<T>(h->h, _char_to_fill(fill_A), n, &one, A, ld_A, x, 1, &zero, y, 1));
            }
            if constexpr(utils::is_complex<T>()) {
                if(op_A == 'N' || op_A == 'H') {
                    T zero = 0.0;
                    T one = 1.0;
                    CHECK(_my_blas_xhemv<T>(h->h, _char_to_fill(fill_A), n, &one, A, ld_A, x, 1, &zero, y, 1));
                }
                else if(op_A == 'C' || op_A == 'T') {
                    utils::remove_complex_t<T> neg_one = -1;
                    if constexpr(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(x) + 1, 2);
                    hemv<T,I>(h, n, A, ld_A, order_A, 'N', fill_A, x, y);
                    if constexpr(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(x) + 1, 2);
                    if constexpr(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(y) + 1, 2);
                }
                else {
                    eslog::error("invalid op_A '%c'\n", op_A);
                }
            }
        }
        else if(order_A == 'R') {
            char fill_A_compl = mgm::fill_change(fill_A);
            char op_A_compl = mgm::operation_combine(op_A, 'C');
            hemv<T,I>(h, n, A, ld_A, 'C', op_A_compl, fill_A_compl, x, y);
        }
        else {
            eslog::error("invalid order_A '%c'\n", order_A);
        }
    }

    template<typename T, typename I>
    void gemv(handle & h, I m, I n, T * A, I ld_A, char order_A, char op_A, T * x, T * y)
    {
        if(order_A == 'C') {
            T zero = 0.0;
            T one = 1.0;
            CHECK(_my_blas_xgemv<T>(h->h, _char_to_operation(op_A), m, n, &one, A, ld_A, x, 1, &zero, y, 1));
        }
        else if(order_A == 'R') {
            char op_A_compl = mgm::operation_combine(op_A, 'T');
            gemv<T,I>(h, m, n, A, ld_A, 'C', op_A_compl, x, y);
        }
        else {
            eslog::error("invalid order_A '%c'\n", order_A);
        }
    }

}
}
}

#include "gpu/gpu_dnblas.inst.hpp"

#endif
#endif
