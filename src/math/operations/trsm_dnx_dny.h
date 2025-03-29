
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
    trsm_dnx_dny() = default;
    trsm_dnx_dny(const trsm_dnx_dny &) = delete;
    trsm_dnx_dny(trsm_dnx_dny &&) = default;
    trsm_dnx_dny & operator=(const trsm_dnx_dny &) = delete;
    trsm_dnx_dny & operator=(trsm_dnx_dny &&) = default;
    ~trsm_dnx_dny() = default;
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
