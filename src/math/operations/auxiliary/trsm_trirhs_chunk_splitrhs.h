
#ifndef SRC_MATH_OPERATIONS_AUXILIARY_TRSM_TRIRHS_CHUNK_SPLITRHS_H
#define SRC_MATH_OPERATIONS_AUXILIARY_TRSM_TRIRHS_CHUNK_SPLITRHS_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_view_new.h"
#include "math/operations/submatrix_csx_csy.h"
#include "math/operations/submatrix_csx_dny.h"
#include "math/operations/submatrix_dnx_dnx_view.h"
#include "math/operations/trsm_csx_dny_dny.h"
#include "math/operations/trsm_dnx_dny.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class trsm_trirhs_chunk_splitrhs
{
public:
    struct config
    {
        char factor_spdn = '_';
        char factor_order = '_';
    };
public:
    trsm_trirhs_chunk_splitrhs() = default;
    trsm_trirhs_chunk_splitrhs(const trsm_trirhs_chunk_splitrhs &) = delete;
    trsm_trirhs_chunk_splitrhs(trsm_trirhs_chunk_splitrhs &&) = default;
    trsm_trirhs_chunk_splitrhs & operator=(const trsm_trirhs_chunk_splitrhs &) = delete;
    trsm_trirhs_chunk_splitrhs & operator=(trsm_trirhs_chunk_splitrhs &&) = default;
    ~trsm_trirhs_chunk_splitrhs() = default;
public:
    void set_config(config cfg_);
    void set_range(size_t rhs_start_, size_t rhs_end_);
    void set_L(MatrixCsxView_new<T,I> * L_);
    void set_X(MatrixDenseView_new<T> * X_);
    void set_X_colpivots(VectorDenseView_new<I> * X_colpivots_);
    void preprocess();
    void perform();
private:
    size_t rhs_start = 0;
    size_t rhs_end = 0;
    size_t rhs_size = 0;
    size_t k_start = 0;
    size_t k_end = 0;
    size_t k_size = 0;
    MatrixCsxView_new<T,I> * L = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
    VectorDenseView_new<I> * X_colpivots = nullptr;
    config cfg;
    bool set_config_called = false;
    bool set_range_called = false;
    bool preprocess_called = false;
    MatrixCsxData_new<T,I> sub_L_sp;
    MatrixDenseData_new<T> sub_L_dn;
    MatrixDenseView_new<T> sub_X;
    submatrix_csx_csy<T,I> op_submatrix_L_sp;
    submatrix_csx_dny<T,I> op_submatrix_L_dn;
    submatrix_dnx_dnx_view<T> op_submatrix_X;
    trsm_csx_dny_dny<T,I> op_trsm_sp;
    trsm_dnx_dny<T> op_trsm_dn;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_AUXILIARY_TRSM_TRIRHS_CHUNK_SPLITRHS_H */
