
#ifndef SRC_MATH_OPERATIONS_COPY_DENSE_H
#define SRC_MATH_OPERATIONS_COPY_DENSE_H

#include "math/primitives_new/matrix_dense_view_new.h"



template<typename T>
class copy_dense
{
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
};



#endif /* SRC_MATH_OPERATIONS_COPY_DENSE_H */
