
#ifndef SRC_MATH_OPERATIONS_CONCAT_CSX_H
#define SRC_MATH_OPERATIONS_CONCAT_CSX_H

#include "math/primitives_new/matrix_csx_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class concat_csx
{
    // ignores output matrix uplo
public:
    concat_csx() = default;
    concat_csx(const concat_csx &) = delete;
    concat_csx(concat_csx &&) = default;
    concat_csx & operator=(const concat_csx &) = delete;
    concat_csx & operator=(concat_csx &&) = default;
    ~concat_csx() = default;
public:
    void set_matrices_src(MatrixCsxView_new<T,I> * A11_, MatrixCsxView_new<T,I> * A12_, MatrixCsxView_new<T,I> * A21_, MatrixCsxView_new<T,I> * A22_);
    void set_matrix_dst(MatrixCsxView_new<T,I> * A_);
    void perform();
    static void do_all(MatrixCsxView_new<T,I> * A11, MatrixCsxView_new<T,I> * A12, MatrixCsxView_new<T,I> * A21, MatrixCsxView_new<T,I> * A22, MatrixCsxView_new<T,I> * A);
private:
    MatrixCsxView_new<T,I> * A11 = nullptr;
    MatrixCsxView_new<T,I> * A12 = nullptr;
    MatrixCsxView_new<T,I> * A21 = nullptr;
    MatrixCsxView_new<T,I> * A22 = nullptr;
    MatrixCsxView_new<T,I> * A = nullptr;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_CONCAT_CSX_H */
