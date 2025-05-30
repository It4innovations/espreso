
#ifndef SRC_MATH_OPERATIONS_COMPLETE_DNX_DNX_H
#define SRC_MATH_OPERATIONS_COMPLETE_DNX_DNX_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class complete_dnx_dnx
{
public:
    complete_dnx_dnx() = default;
    complete_dnx_dnx(const complete_dnx_dnx &) = delete;
    complete_dnx_dnx(complete_dnx_dnx &&) = default;
    complete_dnx_dnx & operator=(const complete_dnx_dnx &) = delete;
    complete_dnx_dnx & operator=(complete_dnx_dnx &&) = default;
    ~complete_dnx_dnx() = default;
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_conj(bool do_conj_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    bool do_conj = true;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_COMPLETE_DNX_DNX_H */
