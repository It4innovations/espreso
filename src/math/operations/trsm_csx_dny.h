
#ifndef SRC_MATH_OPERATIONS_TRSM_CSX_DNY_H
#define SRC_MATH_OPERATIONS_TRSM_CSX_DNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class trsm_csx_dny
{
public:
    trsm_csx_dny() = default;
    trsm_csx_dny(const trsm_csx_dny &) = delete;
    trsm_csx_dny(trsm_csx_dny &&) = delete;
    trsm_csx_dny & operator=(const trsm_csx_dny &) = delete;
    trsm_csx_dny & operator=(trsm_csx_dny &&) = delete;
    ~trsm_csx_dny();
public:
    void set_system_matrix(MatrixCsxView_new<T,I> * A_);
    void set_rhs_sol(MatrixDenseView_new<T> * X_);
    void preprocess();
    void perform();
    void finalize();
private:
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
    MatrixDenseData_new<T> Y;
    spblas::handle_trsm handle;
    bool preprocess_called = false;
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_TRSM_CSX_DNY_H */
