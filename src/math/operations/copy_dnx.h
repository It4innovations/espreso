
#ifndef SRC_MATH_OPERATIONS_COPY_DENSE_H
#define SRC_MATH_OPERATIONS_COPY_DENSE_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class copy_dnx
{
public:
    copy_dnx() = default;
    copy_dnx(const copy_dnx &) = delete;
    copy_dnx(copy_dnx &&) = default;
    copy_dnx & operator=(const copy_dnx &) = delete;
    copy_dnx & operator=(copy_dnx &&) = default;
    ~copy_dnx() = default;
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_conj(bool do_conj_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj = false);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    bool do_conj = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_COPY_DENSE_H */
