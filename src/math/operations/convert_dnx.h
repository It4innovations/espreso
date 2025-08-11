
#ifndef SRC_MATH_OPERATIONS_CONVERT_DNX_H
#define SRC_MATH_OPERATIONS_CONVERT_DNX_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class convert_dnx
{
// respects order
// respects uplo
// does NOT respect diag yet
// does NOT respect conj
// takes symmetry into consideration
// ignores definitness
public:
    convert_dnx() = default;
    convert_dnx(const convert_dnx &) = delete;
    convert_dnx(convert_dnx &&) = default;
    convert_dnx & operator=(const convert_dnx &) = delete;
    convert_dnx & operator=(convert_dnx &&) = default;
    ~convert_dnx() = default;
public:
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst);
private:
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_CONVERT_DNX_H */
