
#ifndef SRC_MATH_OPERATIONS_GEMM_DNX_DNY_DNZ_H
#define SRC_MATH_OPERATIONS_GEMM_DNX_DNY_DNZ_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class gemm_dnx_dny_dnz
{
public:
    gemm_dnx_dny_dnz() = default;
    gemm_dnx_dny_dnz(const gemm_dnx_dny_dnz &) = delete;
    gemm_dnx_dny_dnz(gemm_dnx_dny_dnz &&) = default;
    gemm_dnx_dny_dnz & operator=(const gemm_dnx_dny_dnz &) = delete;
    gemm_dnx_dny_dnz & operator=(gemm_dnx_dny_dnz &&) = default;
    ~gemm_dnx_dny_dnz() = default;
public:
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void set_conj(bool conj_A_, bool conj_B_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta, bool conj_A = false, bool conj_B = false);
private:
    MatrixDenseView_new<T> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    T alpha = T{1};
    T beta = T{0};
    bool conj_A = false;
    bool conj_B = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_GEMM_DNX_DNY_DNZ_H */
