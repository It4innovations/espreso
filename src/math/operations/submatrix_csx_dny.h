
#ifndef SRC_MATH_OPERATIONS_SUBMATRIX_CSX_DNY_H
#define SRC_MATH_OPERATIONS_SUBMATRIX_CSX_DNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class submatrix_csx_dny
{
public:
    void set_matrix_src(const MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(const MatrixDenseView_new<T> * M_dst_);
    void set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_);
    void perform_zerofill();
    void perform_copyvals();
    void perform_all();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    size_t row_start = 0;
    size_t row_end = 0;
    size_t col_start = 0;
    size_t col_end = 0;
    size_t num_rows = 0;
    size_t num_cols = 0;
    bool bound_set = false;
    bool zerofill_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SUBMATRIX_CSX_DNY_H */
