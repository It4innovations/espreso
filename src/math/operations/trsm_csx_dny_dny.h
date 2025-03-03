
#ifndef SRC_MATH_OPERATIONS_TRSM_CSX_DNY_DNY_H
#define SRC_MATH_OPERATIONS_TRSM_CSX_DNY_DNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class trsm_csx_dny_dny
{
// supports in-place (X == B)
public:
    void set_system_matrix(MatrixCsxView_new<T,I> * A_);
    void set_rhs_matrix(MatrixDenseView_new<T> * B_);
    void set_solution_matrix(MatrixDenseView_new<T> * X_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * X);
private:
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_TRSM_CSX_DNY_DNY_H */
