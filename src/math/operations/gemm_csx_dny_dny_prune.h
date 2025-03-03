
#ifndef SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_PRUNE_H
#define SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_PRUNE_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/operations/gemm_csx_dny_dny.h"
#include "math/operations/gemm_dnx_dny_dnz.h"
#include "math/operations/prune_csx_matx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class gemm_csx_dny_dny_prune
{
public:
    gemm_csx_dny_dny_prune() = default;
    gemm_csx_dny_dny_prune(const gemm_csx_dny_dny_prune &) = delete;
    gemm_csx_dny_dny_prune(gemm_csx_dny_dny_prune &&) = delete;
    gemm_csx_dny_dny_prune & operator=(const gemm_csx_dny_dny_prune &) = delete;
    gemm_csx_dny_dny_prune & operator=(gemm_csx_dny_dny_prune &&) = delete;
    ~gemm_csx_dny_dny_prune();
public:
    void set_config(char spdn_A_, bool prune_rows_, bool prune_cols_);
    void set_matrix_A(MatrixCsxView_new<T,I> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void preprocess();
    void perform();
    void finalize();
private:
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    VectorDenseData_new<I> pruned_rows;
    VectorDenseData_new<I> pruned_cols;
    MatrixCsxData_new<T,I> A_pruned_sp;
    MatrixDenseData_new<T> A_pruned_dn;
    MatrixDenseData_new<T> B_pruned;
    MatrixDenseData_new<T> C_pruned;
    MatrixDenseView_new<T> * B_to_use;
    MatrixDenseView_new<T> * C_to_use;
    gemm_csx_dny_dny<T,I> op_gemm_sp;
    gemm_dnx_dny_dnz<T> op_gemm_dn;
    prune_csx_matx<T,I> op_prune_A;
    T alpha = T{1};
    T beta = T{0};
    size_t m = 0;
    size_t n = 0;
    size_t k = 0;
    char spdn_A = '_';
    bool prune_rows = false;
    bool prune_cols = false;
    bool set_config_called = false;
    bool preprocess_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_PRUNE_H */
