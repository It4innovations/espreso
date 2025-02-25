
#ifndef SRC_MATH_OPERATIONS_AUXILIARY_TRSM_TRIRHS_CHUNK_SPLITRHS_H
#define SRC_MATH_OPERATIONS_AUXILIARY_TRSM_TRIRHS_CHUNK_SPLITRHS_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



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
    trsm_trirhs_chunk_splitrhs(trsm_trirhs_chunk_splitrhs &&) = delete;
    trsm_trirhs_chunk_splitrhs & operator=(const trsm_trirhs_chunk_splitrhs &) = delete;
    trsm_trirhs_chunk_splitrhs & operator=(trsm_trirhs_chunk_splitrhs &&) = delete;
    ~trsm_trirhs_chunk_splitrhs();
public:
    void set_config(config cfg_);
    void set_range(size_t rhs_start, size_t rhs_end);
    void set_L(MatrixCsxView_new<T,I> * L_);
    void set_X(MatrixDenseView_new<T> * X_);
    void set_B_colpivots(VectorDenseView_new<T> * B_colpivots_);
    void preprocess();
    void perform();
    void finalize();
private:
    size_t rhs_start = 0;
    size_t rhs_end = 0;
    size_t rhs_size = 0;
    size_t k_start = 0;
    size_t k_end = 0;
    size_t k_size = 0;
    MatrixCsxView_new<T,I> * L = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
    VectorDenseView_new<size_t> * B_colpivots = nullptr;
    config cfg;
    bool set_config_called = false;
    bool set_range_called = false;
    bool preprocess_called = false;
    union {
        MatrixCsxData_new<T,I> sp;
        MatrixDenseData_new<T> dn;
    } sub_L;
    MatrixDenseView_new<T> sub_X;
    submatrix_csx_csy<T,I> op_submatrix_L_sp;
    trsm_csx_dny<T,I> op_trsm_sp;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_AUXILIARY_TRSM_TRIRHS_CHUNK_SPLITRHS_H */
