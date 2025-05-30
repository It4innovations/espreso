
#ifndef SRC_MATH_OPERATIONS_PERMUTE_DNX_DNX_H
#define SRC_MATH_OPERATIONS_PERMUTE_DNX_DNX_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/permutation_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class permute_dnx_dnx
{
    // respects uplo
    // supports in-place (allocates internal tmp matrix)
public:
    permute_dnx_dnx() = default;
    permute_dnx_dnx(const permute_dnx_dnx &) = delete;
    permute_dnx_dnx(permute_dnx_dnx &&) = default;
    permute_dnx_dnx & operator=(const permute_dnx_dnx &) = delete;
    permute_dnx_dnx & operator=(permute_dnx_dnx &&) = default;
    ~permute_dnx_dnx() = default;
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_perm_vector_rows(PermutationView_new<I> * perm_rows_);
    void set_perm_vector_cols(PermutationView_new<I> * perm_cols_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, PermutationView_new<I> * perm_rows, PermutationView_new<I> * perm_cols);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    PermutationView_new<I> * perm_rows = nullptr;
    PermutationView_new<I> * perm_cols = nullptr;
private:
    void perform_primary(PermutationView_new<I> & perm);
    void perform_secdary(PermutationView_new<I> & perm);
    void perform_both(PermutationView_new<I> & perm_primary, PermutationView_new<I> & perm_secdary);
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_PERMUTE_DNX_DNX_H */
