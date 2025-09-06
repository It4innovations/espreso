
#ifndef SRC_MATH_OPERATIONS_SUPERMATIX_DNX_DNX_NONCONTIG_H
#define SRC_MATH_OPERATIONS_SUPERMATIX_DNX_DNX_NONCONTIG_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class supermatrix_dnx_dnx_noncontig
{
    // ignores uplo
public:
    enum struct mode {
        assign,
        accumulate,
        accumulate_atomic,
    };
public:
    supermatrix_dnx_dnx_noncontig() = default;
    supermatrix_dnx_dnx_noncontig(const supermatrix_dnx_dnx_noncontig &) = delete;
    supermatrix_dnx_dnx_noncontig(supermatrix_dnx_dnx_noncontig &&) = default;
    supermatrix_dnx_dnx_noncontig & operator=(const supermatrix_dnx_dnx_noncontig &) = delete;
    supermatrix_dnx_dnx_noncontig & operator=(supermatrix_dnx_dnx_noncontig &&) = default;
    ~supermatrix_dnx_dnx_noncontig() = default;
public:
    void set_matrix_source(MatrixDenseView_new<T> * M_src_);
    void set_matrix_destination(MatrixDenseView_new<T> * M_dst_);
    void set_row_map(VectorDenseView_new<I> * row_map_);
    void set_col_map(VectorDenseView_new<I> * col_map_);
    void set_mode(mode mode_val_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, VectorDenseView_new<I> * row_map, VectorDenseView_new<I> * col_map, mode mode_val = mode::assign);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    VectorDenseView_new<I> * row_map = nullptr;
    VectorDenseView_new<I> * col_map = nullptr;
    mode mode_val = mode::assign;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SUPERMATIX_DNX_DNX_NONCONTIG_H */
