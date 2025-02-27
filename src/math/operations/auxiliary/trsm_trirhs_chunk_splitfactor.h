
#ifndef SRC_MATH_OPERATIONS_AUXILIARY_TRSH_TRIRHS_CHUNK_SPLITFACTOR_H
#define SRC_MATH_OPERATIONS_AUXILIARY_TRSH_TRIRHS_CHUNK_SPLITFACTOR_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_view_new.h"
#include "math/operations/gemm_csx_dny_dny.h"
#include "math/operations/gemm_csx_dny_dny_prune.h"
#include "math/operations/submatrix_csx_csy.h"
#include "math/operations/trsm_csx_dny.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class trsm_trirhs_chunk_splitfactor
{
public:
    struct config
    {
        char trsm_factor_spdn = '_'; // Sparse, Dense
        char trsm_factor_order = '_'; // Rowmajor, Colmajor
        char gemm_factor_spdn = '_'; // Sparse, Dense
        char gemm_factor_order = '_'; // Rowmajor, Colmajor
        char gemm_factor_prune = '_'; // Yes, No
    };
public:
    trsm_trirhs_chunk_splitfactor() = default;
    trsm_trirhs_chunk_splitfactor(const trsm_trirhs_chunk_splitfactor &) = delete;
    trsm_trirhs_chunk_splitfactor(trsm_trirhs_chunk_splitfactor &&) = delete;
    trsm_trirhs_chunk_splitfactor & operator=(const trsm_trirhs_chunk_splitfactor &) = delete;
    trsm_trirhs_chunk_splitfactor & operator=(trsm_trirhs_chunk_splitfactor &&) = delete;
    ~trsm_trirhs_chunk_splitfactor();
public:
    void set_config(config cfg_);
    void set_range(size_t k_start_, size_t k_end_);
    void set_L(MatrixCsxView_new<T,I> * L_);
    void set_X(MatrixDenseView_new<T> * X_);
    void set_X_rowtrails(VectorDenseView_new<I> * X_rowtrails_);
    void preprocess();
    void perform();
    void finalize();
private:
    size_t k_start = 0;
    size_t k_end = 0;
    size_t k_size = 0;
    MatrixCsxView_new<T,I> * L = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
    VectorDenseView_new<I> * X_rowtrails = nullptr;
    config cfg;
    bool set_config_called = false;
    bool set_range_called = false;
    bool preprocess_called = false;
    MatrixCsxData_new<T,I> sub_L_top_sp;
    MatrixDenseData_new<T> sub_L_top_dn;
    MatrixCsxData_new<T,I> sub_L_bot_sp;
    MatrixDenseData_new<T> sub_L_bot_dn;
    MatrixDenseView_new<T> sub_X_top;
    MatrixDenseView_new<T> sub_X_bot;
    submatrix_csx_csy<T,I> op_submatrix_L_top_sp;
    submatrix_csx_csy<T,I> op_submatrix_L_bot_sp;
    trsm_csx_dny<T,I> op_trsm_sp;
    gemm_csx_dny_dny<T,I> op_gemm_normal_sp;
    gemm_csx_dny_dny_prune<T,I> op_gemm_prune;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_AUXILIARY_TRSH_TRIRHS_CHUNK_SPLITFACTOR_H */
