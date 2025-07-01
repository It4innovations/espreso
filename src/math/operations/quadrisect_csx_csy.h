
#ifndef SRC_MATH_OPERATIONS_QUADRISECT_CSX_CSY_H
#define SRC_MATH_OPERATIONS_QUADRISECT_CSX_CSY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/operations/submatrix_csx_csy.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class quadrisect_csx_csy
{
public:
    quadrisect_csx_csy() = default;
    quadrisect_csx_csy(const quadrisect_csx_csy &) = delete;
    quadrisect_csx_csy(quadrisect_csx_csy &&) = default;
    quadrisect_csx_csy & operator=(const quadrisect_csx_csy &) = delete;
    quadrisect_csx_csy & operator=(quadrisect_csx_csy &&) = default;
    ~quadrisect_csx_csy() = default;
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_bounds(size_t row_cut_, size_t col_cut_);
    void setup();
    size_t get_output_matrix_11_nnz();
    size_t get_output_matrix_12_nnz();
    size_t get_output_matrix_21_nnz();
    size_t get_output_matrix_22_nnz();
    void set_matrices_dst(MatrixCsxView_new<T,I> * M_dst_11_, MatrixCsxView_new<T,I> * M_dst_12_, MatrixCsxView_new<T,I> * M_dst_21_, MatrixCsxView_new<T,I> * M_dst_22_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst_11, MatrixCsxView_new<T,I> * M_dst_12, MatrixCsxView_new<T,I> * M_dst_21, MatrixCsxView_new<T,I> * M_dst_22, size_t row_cut, size_t col_cut);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst_11 = nullptr;
    MatrixCsxView_new<T,I> * M_dst_12 = nullptr;
    MatrixCsxView_new<T,I> * M_dst_21 = nullptr;
    MatrixCsxView_new<T,I> * M_dst_22 = nullptr;
    size_t row_cut = 0;
    size_t col_cut = 0;
    bool called_set_bounds = false;
    bool called_setup = false;
    bool called_set_dst = false;
    bool checked = false;
private:
    submatrix_csx_csy<T,I> op_sub_11;
    submatrix_csx_csy<T,I> op_sub_12;
    submatrix_csx_csy<T,I> op_sub_21;
    submatrix_csx_csy<T,I> op_sub_22;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_QUADRISECT_CSX_CSY_H */
