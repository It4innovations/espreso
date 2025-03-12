
#ifndef SRC_MATH_OPERATIONS_PRUNE_CSX_MATX_H
#define SRC_MATH_OPERATIONS_PRUNE_CSX_MATX_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/operations/pruning_subset_csx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class prune_csx_matx
{
public:
    prune_csx_matx() = default;
    prune_csx_matx(const prune_csx_matx &) = delete;
    prune_csx_matx(prune_csx_matx &&) = delete;
    prune_csx_matx & operator=(const prune_csx_matx &) = delete;
    prune_csx_matx & operator=(prune_csx_matx &&) = delete;
    ~prune_csx_matx();
public:
    void set_pruning_mode(bool prune_rows_, bool prune_cols_);
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void setup();
    size_t get_dst_matrix_nrows();
    size_t get_dst_matrix_ncols();
    void set_vector_pruned_rows(VectorDenseView_new<I> * pruned_rows_);
    void set_vector_pruned_cols(VectorDenseView_new<I> * pruned_cols_);
    void preprocess2();
    void set_matrix_dst_sp(MatrixCsxView_new<T,I> * M_dst_sp_);
    void set_matrix_dst_dn(MatrixDenseView_new<T> * M_dst_dn_);
    void perform();
    void finalize();
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst_sp = nullptr;
    MatrixDenseView_new<T> * M_dst_dn = nullptr;
    VectorDenseView_new<I> * pruned_rows_vec;
    VectorDenseView_new<I> * pruned_cols_vec;
    VectorDenseView_new<I> * pruned_idxs_primary = nullptr;
    VectorDenseView_new<I> * pruned_idxs_secdary = nullptr;
    VectorDenseData_new<I> pruned_idxs_secdary_inverse;
    pruning_subset_csx<T,I> op_pruning_subset;
    size_t pruned_nrows = 0;
    size_t pruned_ncols = 0;
    size_t pruned_size_primary = 0;
    size_t pruned_size_secdary = 0;
    bool prune_rows = false;
    bool prune_cols = false;
    bool prune_primary = false;
    bool prune_secdary = false;
    bool called_set_pruning_mode = false;
    bool called_setup = false;
    bool called_preprocess2 = false;
private:
    void perform_sparse();
    void perform_dense();
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_PRUNE_CSX_MATX_H */
