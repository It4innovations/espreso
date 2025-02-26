
#ifndef SRC_MATH_OPERATIONS_FILL_DNX_H
#define SRC_MATH_OPERATIONS_FILL_DNX_H

#include "math/primitivew_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class fill_dnx
{
// respect uplo
public:
    void set_matrix(MatrixDenseView_new<T> * M_);
    void set_value(T val_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M, T val);
private:
    MatrixDenseView_new<T> * M;
    T val = 0;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_FILL_DNX_H */
