
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
    void set_system_matrix(MatrixDenseData_new<T> * A_);
    void set_rhs_sol(MatrixDenseData_new<T> * X_);
    void perform();
    static void do_all(MatrixDenseData_new<T> * A, MatrixDenseData_new<T> * X);
private:
    MatrixDenseData_new<T> * A;
    MatrixDenseData_new<T> * X;
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_TRSM_DNX_DNY_H */
