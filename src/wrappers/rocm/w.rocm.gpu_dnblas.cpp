
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
    void trsm(handle & h, char side, char fill, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X)
    {
        T one = 1.0;
        CHECK(_my_blas_xtrsm(h->h, _char_to_side(side), _char_to_fill(fill), _char_to_operation(transpose), rocblas_diagonal_non_unit, nrows_X, ncols_X, &one, A, ld_A, rhs_sol, ld_X));
    }

    template<typename T, typename I>
    void herk(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C)
    {
        if(utils::is_real<T>() && transpose == 'H') transpose = 'T';
        utils::remove_complex_t<T> zero = 0.0;
        utils::remove_complex_t<T> one = 1.0;
        if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsyrk(h->h, _char_to_fill(out_fill), _char_to_operation(transpose), n, k, &one, A, ld_A, &zero, C, ld_C));
        if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xherk(h->h, _char_to_fill(out_fill), _char_to_operation(transpose), n, k, &one, A, ld_A, &zero, C, ld_C));
    }

    template<typename T, typename I>
    void hemv(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out)
    {
        T zero = 0.0;
        T one = 1.0;
        if constexpr(utils::is_real<T>())    CHECK(_my_blas_xsymv(h->h, _char_to_fill(fill), n, &one, A, ld_A, vec_in, 1, &zero, vec_out, 1));
        if constexpr(utils::is_complex<T>()) CHECK(_my_blas_xhemv(h->h, _char_to_fill(fill), n, &one, A, ld_A, vec_in, 1, &zero, vec_out, 1));
    }

}
}
}

#include "gpu/gpu_dnblas_inst.hpp"

#endif
