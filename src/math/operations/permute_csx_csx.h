
#ifndef SRC_MATH_OPERATIONS_PERMUTE_CSX_CSX_H
#define SRC_MATH_OPERATIONS_PERMUTE_CSX_CSX_H

#include "math/primitivew_new/matrix_csx_view_new.h"
#include "math/primitivew_new/permutation_view_new.h"



template<typename T, typename I>
class permute_csx_csx
{
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void set_perm_vector_rows(PermutationView_new<I> * perm_rows_);
    void set_perm_vector_cols(PermutationView_new<I> * perm_cols_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst, PermutationView_new<I> * perm_rows, PermutationView_new<I> * perm_cols);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    PermutationView_new<I> * perm_rows = nullptr;
    PermutationView_new<I> * perm_cols = nullptr;
private:
    void perform_primary();
    void perform_secdary();
    void perform_both();
};



#endif /* SRC_MATH_OPERATIONS_PERMUTE_CSX_CSX_H */
