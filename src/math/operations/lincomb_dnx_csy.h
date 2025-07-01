
#ifndef SRC_MATH_OPERATIONS_LINCOMB_DNX_CSY_H
#define SRC_MATH_OPERATIONS_LINCOMB_DNX_CSY_H

#include "math/primitives_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class lincomb_dnx_csy
{
// X = alpha * A + beta * B
// inplace allowed (X=A)
// respect uplo
public:
    lincomb_dnx_csy() = default;
    lincomb_dnx_csy(const lincomb_dnx_csy &) = delete;
    lincomb_dnx_csy(lincomb_dnx_csy &&) = default;
    lincomb_dnx_csy & operator=(const lincomb_dnx_csy &) = delete;
    lincomb_dnx_csy & operator=(lincomb_dnx_csy &&) = default;
    ~lincomb_dnx_csy() = default;
public:
    void set_matrix_X(MatrixDenseView_new<T> * X_);
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_B(MatrixCsxView_new<T,I> * B_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * X, T alpha, MatrixDenseView_new<T> * A, T beta, MatrixCsxView_new<T,I> * B);
private:
    MatrixDenseView_new<T> * X = nullptr;
    MatrixDenseView_new<T> * A = nullptr;
    MatrixCsxView_new<T,I> * B = nullptr;
    T alpha = T{0};
    T beta = T{0};
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_LINCOMB_DNX_CSY_H */
