
#ifndef SRC_MATH_OPERATIONS_TRSM_DNX_DNY_H
#define SRC_MATH_OPERATIONS_TRSM_DNX_DNY_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class trsm_dnx_dny
{
public:
    void set_system_matrix(MatrixDenseView_new<T> * A_);
    void set_rhs_sol(MatrixDenseView_new<T> * X_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * X);
private:
    MatrixDenseView_new<T> * A = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_TRSM_DNX_DNY_H */
