
#ifndef SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_PRUNE_H
#define SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_PRUNE_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/operations/gemm_csx_dny_dny.h"



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
    void set_config(char spdn_A_);
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
    VectorDenseData_new<size_t> nonempty;
    VectorDenseData_new<I> A_pruned_sp_ptrs;
    MatrixCsxView_new<T,I> A_pruned_sp;
    MatrixDenseData_new<T> A_pruned_dn;
    MatrixDenseData_new<T> B_pruned;
    MatrixDenseData_new<T> C_pruned;
    MatrixDenseView_new<T> * B_to_use;
    MatrixDenseView_new<T> * C_to_use;
    gemm_csx_dny_dny<T,I> op_gemm_sp;
    T alpha = T{1};
    T beta = T{0};
    char spdn_A = '_';
    bool set_config_called = false;
    bool preprocess_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_GEMM_CSX_DNY_DNY_PRUNE_H */
