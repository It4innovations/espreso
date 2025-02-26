
#ifndef SRC_MATH_OPERATIONS_LINCOMB_MATRIX_DNX_H
#define SRC_MATH_OPERATIONS_LINCOMB_MATRIX_DNX_H

#include "math/primitivew_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class lincomb_matrix_dnx
{
// X = alpha * A + beta * B
// inplace allowed (X=A, X=B)
// respect uplo
public:
    void set_matrix_X(MatrixDenseView_new<T> * X_);
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * X, T alpha, MatrixDenseView_new<T> * A, T beta, MatrixDenseView_new<T> * B);
private:
    MatrixDenseView_new<T> * X = nullptr;
    MatrixDenseView_new<T> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    T alpha = T{0};
    T beta = T{0};
private:
    static void perform_zero(MatrixDenseView_new<T> & X);
    static void perform_one(MatrixDenseView_new<T> & X, T alpha, MatrixDenseView_new<T> & A);
    static void perform_two(MatrixDenseView_new<T> & X, T alpha, MatrixDenseView_new<T> & A, T beta, MatrixDenseView_new<T> & B);
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_LINCOMB_MATRIX_DNX_H */
