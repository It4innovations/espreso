
#ifndef SRC_MATH_OPERATIONS_PERMUTE_CSX_DNY_H
#define SRC_MATH_OPERATIONS_PERMUTE_CSX_DNY_H

#include "math/primitives_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class permute_csx_dny
{
// respects uplo
public:
    permute_csx_dny() = default;
    permute_csx_dny(const permute_csx_dny &) = delete;
    permute_csx_dny(permute_csx_dny &&) = default;
    permute_csx_dny & operator=(const permute_csx_dny &) = delete;
    permute_csx_dny & operator=(permute_csx_dny &&) = default;
    ~permute_csx_dny() = default;
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_perm_rows(PermutationView_new<I> * perm_rows_);
    void set_perm_cols(PermutationView_new<I> * perm_cols_);
    void perform_zerofill();
    void perform_copyvals();
    void perform_all();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, PermutationView_new<I> * perm_rows, PermutationView_new<I> * perm_cols);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    PermutationView_new<I> * perm_rows = nullptr;
    PermutationView_new<I> * perm_cols = nullptr;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_PERMUTE_CSX_DNY_H */
