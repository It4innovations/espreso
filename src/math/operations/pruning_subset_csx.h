
#ifndef MATH_OPERATIONS_PRUNING_SUBSET_CSX_H
#define MATH_OPERATIONS_PRUNING_SUBSET_CSX_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/vector_dense_data_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class pruning_subset_csx
{
public:
    pruning_subset_csx() = default;
    pruning_subset_csx(const pruning_subset_csx &) = delete;
    pruning_subset_csx(pruning_subset_csx &&) = default;
    pruning_subset_csx & operator=(const pruning_subset_csx &) = delete;
    pruning_subset_csx & operator=(pruning_subset_csx &&) = default;
    ~pruning_subset_csx() = default;
public:
    void set_matrix(MatrixCsxView_new<T,I> * M_);
    void set_pruning_mode(bool prune_rows_, bool prune_cols_);
    void setup();
    size_t get_pruned_nrows();
    size_t get_pruned_ncols();
    void set_vector_pruned_rows(VectorDenseView_new<I> * nonempty_rows_);
    void set_vector_pruned_cols(VectorDenseView_new<I> * nonempty_cols_);
    void perform();
private:
    MatrixCsxView_new<T,I> * M = nullptr;
    VectorDenseView_new<I> * nonempty_rows = nullptr;
    VectorDenseView_new<I> * nonempty_cols = nullptr;
    size_t pruned_size_primary = 0;
    size_t pruned_size_secdary = 0;
    VectorDenseData_new<I> nnz_per_secdary;
    bool prune_rows = false;
    bool prune_cols = false;
    bool prune_primary = false;
    bool prune_secdary = false;
    bool called_set_pruning_mode = false;
    bool called_setup = false;
};



}
}
}



#endif /* MATH_OPERATIONS_PRUNING_SUBSET_CSX_H */
