
#ifndef SRC_MATH_OPERATIONS_CONVERT_DNX_DNY_H
#define SRC_MATH_OPERATIONS_CONVERT_DNX_DNY_H

#include "math/primitivew_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class convert_dnx_dny
{
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_conj(bool do_conj_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj);
private:
    MatrixDenseView_new<T> * M_src;
    MatrixDenseView_new<T> * M_dst;
    bool do_conj = true;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_CONVERT_DNX_DNY_H */
