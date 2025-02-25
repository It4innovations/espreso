
#ifndef SRC_MATH_OPERATIONS_SUBMATRIX_CSX_CSY_MAP_H
#define SRC_MATH_OPERATIONS_SUBMATRIX_CSX_CSY_MAP_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/vector_dense_data_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class submatrix_csx_csy_map_map
{
public:
    submatrix_csx_csy_map() = default;
    submatrix_csx_csy_map(const submatrix_csx_csy_map &) = delete;
    submatrix_csx_csy_map(submatrix_csx_csy_map &&) = delete;
    submatrix_csx_csy_map & operator=(const submatrix_csx_csy_map &) = delete;
    submatrix_csx_csy_map & operator=(submatrix_csx_csy_map &&) = delete;
    ~submatrix_csx_csy_map();
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_);
    void setup();
    void get_output_matrix_nnz();
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void perform_pattern();
    void perform_values();
    void perform_all();
    void finalize();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, size_t row_start, size_t row_end, size_t col_start, size_t col_end);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
    VectorDenseData_new<I> map;
    size_t row_start = 0;
    size_t row_end = 0;
    size_t col_start = 0;
    size_t col_end = 0;
    size_t num_rows = 0;
    size_t num_cols = 0;
    size_t nnz_output = 0;
    bool set_bounds_called = false;
    bool setup_called = false;
    bool perform_pattern_called = false;
private:
    void perform_pattern_same_order();
    void perform_pattern_diff_order();
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SUBMATRIX_CSX_CSY_MAP_H */
