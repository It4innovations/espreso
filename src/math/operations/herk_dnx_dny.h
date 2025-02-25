
#ifndef SRC_MATH_OPERATIONS_SYRK_DNX_DNY_H
#define SRC_MATH_OPERATIONS_SYRK_DNX_DNY_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class herk_dnx_dny
{
public:
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_mode(blas::herk_mode mode_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * C, herk_mode mode, T alpha, T beta);
private:
    MatrixDenseView_new<T> * A = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    T alpha = T{1};
    T beta = T{0};
    blas::herk_mode mode;
    bool mode_set = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SYRK_DNX_DNY_H */
