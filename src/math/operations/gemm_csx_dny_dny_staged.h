
#ifndef SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_STAGED_H
#define SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_STAGED_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/wrappers/math.spblas.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class gemm_csx_dny_dny_staged
{
public:
    gemm_csx_dny_dny_staged() = default;
    gemm_csx_dny_dny_staged(const gemm_csx_dny_dny_staged &) = delete;
    gemm_csx_dny_dny_staged(gemm_csx_dny_dny_staged &&) = delete;
    gemm_csx_dny_dny_staged & operator=(const gemm_csx_dny_dny_staged &) = delete;
    gemm_csx_dny_dny_staged & operator=(gemm_csx_dny_dny_staged &&) = delete;
    ~gemm_csx_dny_dny_staged();
public:
    void set_matrix_A(MatrixCsxView_new<T,I> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void preprocess();
    void perform();
    void finalize();
private:
    spblas::handle_mm handle_abc;
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    T alpha = T{1};
    T beta = T{0};
    bool preprocess_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_STAGED_H */
