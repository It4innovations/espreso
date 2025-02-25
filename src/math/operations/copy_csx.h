
#ifndef SRC_MATH_OPERATIONS_COPY_CSX_H
#define SRC_MATH_OPERATIONS_COPY_CSX_H

#include "math/primitives_new/matrix_csx_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class copy_csx
{
public:
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst);
private:
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixCsxView_new<T,I> * M_dst = nullptr;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_COPY_CSX_H */
