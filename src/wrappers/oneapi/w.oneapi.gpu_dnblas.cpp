
#ifdef HAVE_ONEAPI

#include "gpu/gpu_dnblas.h"
#include "w.oneapi.gpu_management.h"
#include <oneapi/mkl.hpp>

namespace espreso {
namespace gpu {
namespace dnblas {

    namespace onemkl = oneapi::mkl;
    namespace oneblas = onemkl::blas;

    namespace
    {
        static inline onemkl::uplo _char_to_uplofill(char c)
        {
            switch(c) {
                case 'U': return onemkl::uplo::upper;
                case 'L': return onemkl::uplo::lower;
                default: eslog::error("invalid fill '%c'\n", c);
            }
        }

        static inline onemkl::transpose _char_to_operation(char c)
        {
            switch(c) {
                case 'N': return onemkl::transpose::nontrans;
                case 'T': return onemkl::transpose::trans;
                case 'H': return onemkl::transpose::conjtrans;
                default: eslog::error("invalid operation '%c'\n", c);
            }
        }

        static inline onemkl::side _char_to_side(char c)
        {
            switch(c) {
                case 'L': return onemkl::side::left;
                case 'R': return onemkl::side::right;
                default: eslog::error("invalid side '%c'\n", c);
            }
        }
    }

    struct _handle
    {
        mgm::queue q;
        sycl::queue qq;
    };

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        h->q = q;
        h->qq = q->q;
    }

    void handle_destroy(handle & h)
    {
        h.reset();
    }

    void buffer_collect_size(handle & h, size_t & buffersize, const std::function<void(void)> & f)
    {
        f();
        buffersize = 0;
    }

    void buffer_set(handle & h, void * ptr, size_t size) {}

    void buffer_unset(handle & h) {}

    template<typename T, typename I>
    void trsv(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x)
    {
        if(op_A == 'C') {
            char order_A_compl = mgm::order_change(order_A);
            char op_A_compl = mgm::operation_combine(op_A, 'T');
            char fill_A_compl = mgm::fill_change(fill_A);
            trsv<T,I>(h, n, A, ld_A, order_A_compl, op_A_compl, fill_A_compl, x);
        }
        else {
            if(order_A == 'R') {
                oneblas::row_major::trsv(h->qq, _char_to_uplofill(fill_A), _char_to_operation(op_A), onemkl::diag::nonunit, n, A, ld_A, x, 1);
            }
            else if(order_A == 'C') {
                oneblas::column_major::trsv(h->qq, _char_to_uplofill(fill_A), _char_to_operation(op_A), onemkl::diag::nonunit, n, A, ld_A, x, 1);
            }
            else {
                eslog::error("invalid order_A '%c'\n", order_A);
            }
        }
    }

    template<typename T, typename I>
    void trsm(handle & h, char side, I n, I nrhs, T * A, I ld_A, char order_A, char op_A, char fill_A, T * X, I ld_X, char order_X, char op_X)
    {
        if(order_X == order_A) {
            if(op_X == 'N') {
                if(op_A == 'C') {
                    char order_A_compl = mgm::order_change(order_A);
                    char op_A_compl = mgm::operation_combine(op_A, 'T');
                    char fill_A_compl = mgm::fill_change(fill_A);
                    trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A_compl, op_A_compl, fill_A_compl, X, ld_X, order_X, op_X);
                }
                else {
                    T one = 1.0;
                    I nrows_X = (side == 'L' ? n : nrhs);
                    I ncols_X = (side == 'L' ? nrhs : n);
                    if(order_A == 'R') {
                        oneblas::row_major::trsm(h->qq, _char_to_side(side), _char_to_uplofill(fill_A), _char_to_operation(op_A), onemkl::diag::nonunit, nrows_X, ncols_X, one, A, ld_A, X, ld_X);
                    }
                    else if(order_A == 'C') {
                        oneblas::column_major::trsm(h->qq, _char_to_side(side), _char_to_uplofill(fill_A), _char_to_operation(op_A), onemkl::diag::nonunit, nrows_X, ncols_X, one, A, ld_A, X, ld_X);
                    }
                    else {
                        eslog::error("invalid order_A '%c'\n", order_A);
                    }
                }
            }
            else if(op_X == 'C') {
                char op_X_compl = 'N';
                char op_A_compl = mgm::operation_combine(op_A, 'C');
                trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, op_A_compl, fill_A, X, ld_X, order_X, op_X_compl);
            }
            else if(op_X == 'T' || op_X == 'H') {
                char op_X_compl = 'N';
                char op_A_compl = mgm::operation_combine(op_A, op_X);
                char side_compl = mgm::side_change(side);
                trsm<T,I>(h, side_compl, n, nrhs, A, ld_A, order_A, op_A_compl, fill_A, X, ld_X, order_X, op_X_compl);
            }
            else {
                eslog::error("invalid op_X %c\n", op_X);
            }
        }
        else {
            char order_X_compl = order_A;
            char op_X_compl = mgm::operation_combine(op_X, 'T');
            trsm<T,I>(h, side, n, nrhs, A, ld_A, order_A, op_A, fill_A, X, ld_X, order_X_compl, op_X_compl);
        }
    }

    template<typename T, typename I>
    void herk(handle & h, I n, I k, T * A, I ld_A, char order_A, char op_A, T * C, I ld_C, char order_C, char fill_C)
    {
        if(order_C == order_A) {
            T zero = 0.0;
            T one = 1.0;
            if constexpr(utils::is_real<T>()) {
                op_A = mgm::operation_remove_conj(op_A);
            }
            if(order_A == 'R') {
                if constexpr(utils::is_real<T>())    oneblas::row_major::syrk(h->qq, _char_to_uplofill(fill_C), _char_to_operation(op_A), n, k, one, A, ld_A, zero, C, ld_C);
                if constexpr(utils::is_complex<T>()) oneblas::row_major::herk(h->qq, _char_to_uplofill(fill_C), _char_to_operation(op_A), n, k, one, A, ld_A, zero, C, ld_C);
            }
            else if(order_A == 'C') {
                if constexpr(utils::is_real<T>())    oneblas::column_major::syrk(h->qq, _char_to_uplofill(fill_C), _char_to_operation(op_A), n, k, one, A, ld_A, zero, C, ld_C);
                if constexpr(utils::is_complex<T>()) oneblas::column_major::herk(h->qq, _char_to_uplofill(fill_C), _char_to_operation(op_A), n, k, one, A, ld_A, zero, C, ld_C);
            }
            else {
                eslog::error("invalid order_A '%c'\n", order_A);
            }
        }
        else {
            char order_C_compl = order_A;
            char fill_C_compl = mgm::fill_change(fill_C);
            char op_A_compl = mgm::operation_combine(op_A, 'C');
            herk<T,I>(h, n, k, A, ld_A, order_A, op_A_compl, C, ld_C, order_C_compl, fill_C_compl);
        }
    }

    template<typename T, typename I>
    void hemv(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x, T * y)
    {
        if(op_A == 'N') {
            T zero = 0.0;
            T one = 1.0;
            if(order_A == 'R') {
                if constexpr(utils::is_real<T>())    oneblas::row_major::symv(h->qq, _char_to_uplofill(fill_A), n, one, A, ld_A, x, 1, zero, y, 1);
                if constexpr(utils::is_complex<T>()) oneblas::row_major::hemv(h->qq, _char_to_uplofill(fill_A), n, one, A, ld_A, x, 1, zero, y, 1);
            }
            else if(order_A == 'C') {
                if constexpr(utils::is_real<T>())    oneblas::column_major::symv(h->qq, _char_to_uplofill(fill_A), n, one, A, ld_A, x, 1, zero, y, 1);
                if constexpr(utils::is_complex<T>()) oneblas::column_major::hemv(h->qq, _char_to_uplofill(fill_A), n, one, A, ld_A, x, 1, zero, y, 1);
            }
            else {
                eslog::error("invalid order_A '%c'\n", order_A);
            }
        }
        else {
            eslog::error("hemv can't work with op_A '%c'\n", op_A);
        }
    }

    template<typename T, typename I>
    void gemv(handle & h, I m, I n, T * A, I ld_A, char order_A, char op_A, T * x, T * y)
    {
        T zero = 0.0;
        T one = 1.0;
        if(order_A == 'R') {
            oneblas::row_major::gemv(h->qq, _char_to_operation(op_A), m, n, one, A, ld_A, x, 1, zero, y, 1);
        }
        else if(order_A == 'C') {
            oneblas::column_major::gemv(h->qq, _char_to_operation(op_A), m, n, one, A, ld_A, x, 1, zero, y, 1);
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
