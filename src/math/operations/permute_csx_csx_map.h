
#ifndef SRC_MATH_OPERATIONS_PERMUTE_CSX_CSX_MAP_H
#define SRC_MATH_OPERATIONS_PERMUTE_CSX_CSX_MAP_H

#include "math/primitives_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class permute_csx_csx_map
{
// supports uplo
public:
    permute_csx_csx_map() = default;
    permute_csx_csx_map(const permute_csx_csx_map &) = delete;
    permute_csx_csx_map(permute_csx_csx_map &&) = default;
    permute_csx_csx_map & operator=(const permute_csx_csx_map &) = delete;
    permute_csx_csx_map & operator=(permute_csx_csx_map &&) = default;
    ~permute_csx_csx_map() = default;
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void set_perm_rows(PermutationView_new<I> * perm_rows_);
    void set_perm_cols(PermutationView_new<I> * perm_cols_);
    void perform_pattern();
    void perform_values();
    void perform_all();
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    PermutationView_new<I> * perm_rows = nullptr;
    PermutationView_new<I> * perm_cols = nullptr;
    VectorDenseData_new<I> map_dst_to_src;
    bool called_perform_pattern = false;
private:
    void perform_primary(PermutationView_new<I> & perm);
    void perform_secdary(PermutationView_new<I> & perm);
    void perform_both(PermutationView_new<I> & perm_primary, PermutationView_new<I> & perm_secdary);
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_PERMUTE_CSX_CSX_MAP_H */
