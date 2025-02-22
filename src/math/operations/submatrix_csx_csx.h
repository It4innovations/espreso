
#ifndef SRC_MATH_OPERATIONS_SUBMATRIX_CSX_CSX_H
#define SRC_MATH_OPERATIONS_SUBMATRIX_CSX_CSX_H

#include "math/primitives_new/matrix_csx_view_new.h"



template<typename T, typename I>
class submatrix_csx_csx
{
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_);
    void setup();
    void get_output_matrix_nnz();
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    size_t row_start = 0;
    size_t row_end = 0;
    size_t col_start = 0;
    size_t col_end = 0;
    size_t num_rows = 0;
    size_t num_cols = 0;
    size_t nnz_output = 0;
    bool bounds_set = false;
    bool setup_called = false;
};



#endif /* SRC_MATH_OPERATIONS_SUBMATRIX_CSX_CSX_H */
