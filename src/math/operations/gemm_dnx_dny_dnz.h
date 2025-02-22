
#ifndef SRC_MATH_OPERATIONS_GEMM_DNX_DNY_DNZ_H
#define SRC_MATH_OPERATIONS_GEMM_DNX_DNY_DNZ_H

#include "math/primitives_new/matrix_dense_view_new.h"



template<typename T>
class gemm_dnx_dny_dnz
{
public:
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta);
private:
    MatrixDenseView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T,I> * B = nullptr;
    MatrixDenseView_new<T,I> * C = nullptr;
    T alpha = T{1};
    T beta = T{0};
};



#endif /* SRC_MATH_OPERATIONS_GEMM_DNX_DNY_DNZ_H */
