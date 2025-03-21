
#ifndef SRC_MATH_OPERATIONS_SUBMATRIX_DNX_DNX_VIEW_H
#define SRC_MATH_OPERATIONS_SUBMATRIX_DNX_DNX_VIEW_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class submatrix_dnx_dnx_view
{
public:
    submatrix_dnx_dnx_view() = default;
    submatrix_dnx_dnx_view(const submatrix_dnx_dnx_view &) = delete;
    submatrix_dnx_dnx_view(submatrix_dnx_dnx_view &&) = default;
    submatrix_dnx_dnx_view & operator=(const submatrix_dnx_dnx_view &) = delete;
    submatrix_dnx_dnx_view & operator=(submatrix_dnx_dnx_view &&) = default;
    ~submatrix_dnx_dnx_view() = default;
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    size_t row_start = 0;
    size_t row_end = 0;
    size_t col_start = 0;
    size_t col_end = 0;
    size_t num_rows = 0;
    size_t num_cols = 0;
    bool called_set_bounds = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SUBMATRIX_DNX_DNX_VIEW_H */
