
#ifdef HAVE_ROCM

#include "gpu/gpu_dnblas.h"
#include "w.rocm.gpu_management.h"
#include "basis/utilities/utils.h"

#include <rocblas.h>



inline void _check(rocblas_status status, const char *file, int line)
{
    if (status != rocblas_status_success && status != rocblas_status_size_unchanged && status != rocblas_status_size_increased)
    {
        espreso::eslog::error("ROCBLAS Error %d. In file '%s' on line %d\n", status, file, line);
    }
}



namespace espreso {
namespace gpu {
namespace dnblas {

    namespace
    {
        template<typename T> struct cpp_to_rocblas_type { using type = T; };
        template<> struct cpp_to_rocblas_type<std::complex<float>> { using type = rocblas_float_complex ; };
        template<> struct cpp_to_rocblas_type<std::complex<double>> { using type = rocblas_double_complex ; };
        template<typename T> using cpp_to_rocblas_type_t = typename cpp_to_rocblas_type<T>::type;

        template<typename T>
        static rocblas_status _my_blas_xtrsv(rocblas_handle handle, rocblas_fill uplo, rocblas_operation transA, rocblas_diagonal diag, rocblas_int m, const T *A, rocblas_int lda, T *x, rocblas_int incx)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return rocblas_strsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
            if constexpr(std::is_same_v<T,double>)               return rocblas_dtrsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return rocblas_ctrsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<double>>) return rocblas_ztrsv(handle, uplo, transA, diag, m, (U*)A, lda, (U*)x, incx);
        }

        template<typename T>
        static rocblas_status _my_blas_xtrsm(rocblas_handle handle, rocblas_side side, rocblas_fill uplo, rocblas_operation transA, rocblas_diagonal diag, rocblas_int m, rocblas_int n, const T *alpha, const T *A, rocblas_int lda, T *B, rocblas_int ldb)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return rocblas_strsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
            if constexpr(std::is_same_v<T,double>)               return rocblas_dtrsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return rocblas_ctrsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
            if constexpr(std::is_same_v<T,std::complex<double>>) return rocblas_ztrsm(handle, side, uplo, transA, diag, m, n, (U*)alpha, (U*)A, lda, (U*)B, ldb);
        }

        template<typename T>
        static rocblas_status _my_blas_xsyrk(rocblas_handle handle, rocblas_fill uplo, rocblas_operation transA, rocblas_int n, rocblas_int k, const T *alpha, const T *A, rocblas_int lda, const T *beta, T *C, rocblas_int ldc)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,float>)  return rocblas_ssyrk(handle, uplo, transA, n, k, (U*)alpha, (U*)A, lda, (U*)beta, (U*)C, ldc);
            if constexpr(std::is_same_v<T,double>) return rocblas_dsyrk(handle, uplo, transA, n, k, (U*)alpha, (U*)A, lda, (U*)beta, (U*)C, ldc);
        }

        template<typename T>
        static rocblas_status _my_blas_xherk(rocblas_handle handle, rocblas_fill uplo, rocblas_operation transA, rocblas_int n, rocblas_int k, const utils::remove_complex_t<T> *alpha, const T *A, rocblas_int lda, const utils::remove_complex_t<T> *beta, T *C, rocblas_int ldc)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,std::complex<float>>)  return rocblas_cherk(handle, uplo, transA, n, k, alpha, (U*)A, lda, beta, (U*)C, ldc);
            if constexpr(std::is_same_v<T,std::complex<double>>) return rocblas_zherk(handle, uplo, transA, n, k, alpha, (U*)A, lda, beta, (U*)C, ldc);
        }

        template<typename T>
        static rocblas_status _my_blas_xsymv(rocblas_handle handle, rocblas_fill uplo, rocblas_int n, const T *alpha, const T *A, rocblas_int lda, const T *x, rocblas_int incx, const T *beta, T *y, rocblas_int incy)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,float>)  return rocblas_ssymv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,double>) return rocblas_dsymv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
        }

        template<typename T>
        static rocblas_status _my_blas_xhemv(rocblas_handle handle, rocblas_fill uplo, rocblas_int n, const T *alpha, const T *A, rocblas_int lda, const T *x, rocblas_int incx, const T *beta, T *y, rocblas_int incy)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,std::complex<float>>)  return rocblas_csymv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,std::complex<double>>) return rocblas_zsymv(handle, uplo, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
        }

        template<typename T>
        static rocblas_status _my_blas_xgemv(rocblas_handle handle, rocblas_operation trans, int m, int n, const T *alpha, const T *A, int lda, const T *x, int incx, const T *beta, T *y, int incy)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return rocblas_sgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,double>)               return rocblas_dgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return rocblas_cgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
            if constexpr(std::is_same_v<T,std::complex<double>>) return rocblas_zgemv(handle, trans, m, n, (U*)alpha, (U*)A, lda, (U*)x, incx, (U*)beta, (U*)y, incy);
        }

        template<typename T>
        static rocblas_status _my_blas_xscal(rocblas_handle handle, int n, const T *alpha, T *x, int incx)
        {
            using U = cpp_to_rocblas_type_t<T>;
            if constexpr(std::is_same_v<T,float>)                return rocblas_sscal(handle, n, (U*)alpha, (U*)x, incx);
            if constexpr(std::is_same_v<T,double>)               return rocblas_dscal(handle, n, (U*)alpha, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<float>>)  return rocblas_cscal(handle, n, (U*)alpha, (U*)x, incx);
            if constexpr(std::is_same_v<T,std::complex<double>>) return rocblas_zscal(handle, n, (U*)alpha, (U*)x, incx);
        }

        static rocblas_fill _char_to_fill(char c)
        {
            switch(c)
            {
                case 'U': return rocblas_fill_upper;
                case 'L': return rocblas_fill_lower;
                default: eslog::error("invalid fill '%c'\n", c);
            }
        }

        static rocblas_operation _char_to_operation(char c)
        {
            switch(c)
            {
                case 'N': return rocblas_operation_none;
                case 'T': return rocblas_operation_transpose;
                case 'H': return rocblas_operation_conjugate_transpose;
                default: eslog::error("invalid operation '%c'\n", c);
            }
        }

        static rocblas_side _char_to_side(char c)
        {
            switch(c)
            {
                case 'L': return rocblas_side_left;
                case 'R': return rocblas_side_right;
                default: eslog::error("invalid side '%c'\n", c);
            }
        }
    }

    struct _handle
    {
        rocblas_handle h;
    };

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        CHECK(rocblas_create_handle(&h->h));
        CHECK(rocblas_set_stream(h->h, q->stream));
    }

    void handle_destroy(handle & h)
    {
        if(h.get() == nullptr) return;

        CHECK(rocblas_destroy_handle(h->h));
        h.reset();
    }

    void buffer_collect_size(handle & h, size_t & buffersize, const std::function<void(void)> & f)
    {
        CHECK(rocblas_start_device_memory_size_query(h->h));
        f();
        CHECK(rocblas_stop_device_memory_size_query(h->h, &buffersize));
    }

    void buffer_set(handle & h, void * ptr, size_t size)
    {
        CHECK(rocblas_set_workspace(h->h, ptr, size));
    }

    void buffer_unset(handle & h)
    {
        CHECK(rocblas_set_workspace(h->h, nullptr, 0));
    }

    template<typename T, typename I>
    void trsv(handle & h, char fill, char transpose, I n, I ld, T * matrix, T * rhs_sol)
    {
        CHECK(_my_blas_xtrsv(h->h, _char_to_fill(fill), _char_to_operation(transpose), rocblas_diagonal_non_unit, n, matrix, ld, rhs_sol, 1));
    }

    template<typename T, typename I>
    void trsm(handle & h, char side, I n, I nrhs, T * A, I ld_A, char order_A, char op_A, char fill_A, T * X, I ld_X, char order_X, char op_X)
    {
        // eslog::info("TRSM side %c order_A %c op_A %c fill_A %c order_X %c op_X %c\n", side, order_A, op_A, fill_A, order_X, op_X);
        // solve op(A) op(X) = op(B) or op(X) op(A) = op(B), where input B is stored in X
        if(order_A == 'C')
        {
            if(order_X == 'C')
            {
                if(op_X == 'N')
                {
                    if(op_A == 'C') // conjugate_only operation not supported ...
                    {
                        utils::remove_complex_t<T> neg_one = -1;
                        I nvals_X = (side == 'L' ? nrhs : n) * ld_X; // i know order_X == 'C'
                        if(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, nvals_X, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(X) + 1, 2);
                        trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, 'N', fill_A, X, ld_X, order_X, op_X);
                        if(utils::is_complex<T>()) _my_blas_xscal<utils::remove_complex_t<T>>(h->h, nvals_X, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(X) + 1, 2);
                    }
                    else
                    {
                        T one = 1.0;
                        I nrows_X = (side == 'L' ? n : nrhs);
                        I ncols_X = (side == 'L' ? nrhs : n);
                        #pragma omp critical(espreso_gpu_dnblas_rocm)
                        CHECK(_my_blas_xtrsm<T>(h->h, _char_to_side(side), _char_to_fill(fill_A), _char_to_operation(op_A), rocblas_diagonal_non_unit, nrows_X, ncols_X, &one, A, ld_A, X, ld_X));
                    }
                }
                else if(op_X == 'C')
                {
                    char op_A_comb = mgm::operation_combine(op_A, 'C');
                    trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, op_A_comb, fill_A, X, ld_X, order_X, 'N');
                }
                else if(op_X == 'T' || op_X == 'H')
                {
                    char side_compl = mgm::side_change(side);
                    char op_A_comb = mgm::operation_combine(op_A, op_X);
                    trsm<T,I>(h, side_compl, n, nrhs, A, ld_A, order_A, op_A_comb, fill_A, X, ld_X, order_X, 'N');
                }
                else eslog::error("invalid op_X %c\n", op_X);
            }
            else if(order_X == 'R')
            {
                char op_X_compl = mgm::operation_combine(op_X, 'T');
                trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, op_A, fill_A, X, ld_X, 'C', op_X_compl);
            }
            else eslog::error("invalid order_X '%c'\n", order_X);
        }
        else if(order_A == 'R')
        {
            char op_A_compl = mgm::operation_combine(op_A, 'T');
            char fill_A_compl = mgm::fill_change(fill_A);
            trsm<T,I>(h, side, n, nrhs, A, ld_A, 'C', op_A_compl, fill_A_compl, X, ld_X, order_X, op_X);
        }
        else eslog::error("invalid order_A '%c'\n", order_A);
    }

    template<typename T, typename I>
    void herk(handle & h, I n, I k, T * A, I ld_A, char order_A, char op_A, T * C, I ld_C, char order_C, char fill_C)
    {
        // C = op(A) * op(A)^H
        if(order_C == 'C')
        {
            if(order_A == 'C')
            {
                if(utils::is_real<T>() && op_A == 'H') op_A = 'T';
                utils::remove_complex_t<T> zero = 0.0;
                utils::remove_complex_t<T> one = 1.0;
                #pragma omp critical(espreso_gpu_dnblas_rocm)
                if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsyrk<T>(h->h, _char_to_fill(fill_C), _char_to_operation(op_A), n, k, &one, A, ld_A, &zero, C, ld_C));
                #pragma omp critical(espreso_gpu_dnblas_rocm)
                if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xherk<T>(h->h, _char_to_fill(fill_C), _char_to_operation(op_A), n, k, &one, A, ld_A, &zero, C, ld_C));
            }
            else if(order_A == 'R')
            {
                char op_A_compl = mgm::operation_combine(op_A, 'T');
                herk<T,I>(h, n, k, A, ld_A, 'C', op_A_compl, C, ld_C, order_C, fill_C);
            }
            else eslog::error("invalid order_A '%c'\n", order_A);
        }
        else if(order_C == 'R')
        {
            char fill_C_compl = mgm::fill_change(fill_C);
            char op_A_compl = mgm::operation_combine(op_A, 'C');
            herk<T,I>(h, n, k, A, ld_A, order_A, op_A_compl, C, ld_C, 'C', fill_C_compl);
        }
        else eslog::error("invalid order_C '%c'\n", order_C);
    }

    template<typename T, typename I>
    void hemv(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x, T * y)
    {
        if(order_A == 'C')
        {
            if constexpr(utils::is_real<T>())
            {
                T zero = 0.0;
                T one = 1.0;
                CHECK(_my_blas_xsymv<T>(h->h, _char_to_fill(fill_A), n, &one, A, ld_A, x, 1, &zero, y, 1));
            }
            if constexpr(utils::is_complex<T>())
            {
                if(utils::is_real<T>() || op_A == 'N' || op_A == 'H')
                {
                    T zero = 0.0;
                    T one = 1.0;
                    if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xhemv<T>(h->h, _char_to_fill(fill_A), n, &one, A, ld_A, x, 1, &zero, y, 1));
                }
                else if(op_A == 'C' || op_A == 'T')
                {
                    utils::remove_complex_t<T> neg_one = -1;
                    _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(x) + 1, 2);
                    hemv<T,I>(h, n, A, ld_A, order_A, 'N', fill_A, x, y);
                    _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(x) + 1, 2);
                    _my_blas_xscal<utils::remove_complex_t<T>>(h->h, n, &neg_one, reinterpret_cast<utils::remove_complex_t<T>*>(y) + 1, 2);
                }
            }
        }
        else if(order_A == 'R')
        {
            char fill_A_compl = mgm::fill_change(fill_A);
            char op_A_compl = mgm::operation_combine(op_A, 'C');
            hemv<T,I>(h, n, A, ld_A, 'C', op_A_compl, fill_A_compl, x, y);
        }
        else eslog::error("invalid order_A '%c'\n", order_A);
    }

    template<typename T, typename I>
    void gemv(handle & h, I m, I n, T * A, I ld_A, char order_A, char op_A, T * x, T * y)
    {
        if(order_A == 'C')
        {
            T zero = 0.0;
            T one = 1.0;
            CHECK(_my_blas_xgemv<T>(h->h, _char_to_operation(op_A), m, n, &one, A, ld_A, x, 1, &zero, y, 1));
        }
        else if(order_A == 'R')
        {
            char op_A_compl = mgm::operation_combine(op_A, 'T');
            gemv<T,I>(h, m, n, A, ld_A, 'C', op_A_compl, x, y);
        }
        else eslog::error("invalid order_A '%c'\n", order_A);
    }

}
}
}

#include "gpu/gpu_dnblas.inst.hpp"

#endif
