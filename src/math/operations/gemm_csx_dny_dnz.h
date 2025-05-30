
#ifndef SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNZ_H
#define SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNZ_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/wrappers/math.spblas.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class gemm_csx_dny_dnz
{
public:
    gemm_csx_dny_dnz() = default;
    gemm_csx_dny_dnz(const gemm_csx_dny_dnz &) = delete;
    gemm_csx_dny_dnz(gemm_csx_dny_dnz &&) = default;
    gemm_csx_dny_dnz & operator=(const gemm_csx_dny_dnz &) = delete;
    gemm_csx_dny_dnz & operator=(gemm_csx_dny_dnz &&) = default;
    ~gemm_csx_dny_dnz() = default;
public:
    void set_matrix_A(MatrixCsxView_new<T,I> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta);
private:
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    T alpha = T{1};
    T beta = T{0};
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNZ_H */
