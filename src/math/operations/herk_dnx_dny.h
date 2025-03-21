
#ifndef SRC_MATH_OPERATIONS_SYRK_DNX_DNY_H
#define SRC_MATH_OPERATIONS_SYRK_DNX_DNY_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/wrappers/math.blas.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class herk_dnx_dny
{
public:
    using Treal = utils::remove_complex_t<T>;
public:
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_mode(blas::herk_mode mode_);
    void set_coefficients(Treal alpha_, Treal beta_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * C, blas::herk_mode mode, Treal alpha, Treal beta);
private:
    MatrixDenseView_new<T> * A = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    Treal alpha = Treal{1};
    Treal beta = Treal{0};
    blas::herk_mode mode;
    bool mode_set = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SYRK_DNX_DNY_H */
