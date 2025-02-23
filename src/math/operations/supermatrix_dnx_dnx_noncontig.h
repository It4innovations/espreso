
#ifndef SRC_MATH_OPERATIONS_SUPERMATIX_DNX_DNX_NONCONTIG_H
#define SRC_MATH_OPERATIONS_SUPERMATIX_DNX_DNX_NONCONTIG_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



template<typename T>
class supermatrix_dnx_dnx_noncontig
{
public:
    void set_matrix_source(MatrixDenseView_new<T> * M_src_);
    void set_matrix_destinatino(MatrixDenseView_new<T> * M_dst_);
    void set_row_map(VectorDenseView_new<size_t> * row_map_);
    void set_col_map(VectorDenseView_new<size_t> * col_map_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, VectorDenseView_new<size_t> * row_map, VectorDenseView_new<size_t> * col_map);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    VectorDenseView_new<size_t> * row_map = nullptr;
    VectorDenseView_new<size_t> * col_map = nullptr;
};



#endif /* SRC_MATH_OPERATIONS_SUPERMATIX_DNX_DNX_NONCONTIG_H */
