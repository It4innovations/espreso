
#ifndef SRC_MATH_OPERATIONS_COMPLETE_DNX_H
#define SRC_MATH_OPERATIONS_COMPLETE_DNX_H

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
class complete_dnx
{
public:
    void set_matrix(MatrixDenseView_new<T> * M_);
    void set_orig_uplo(char orig_uplo_);
    void set_conj(bool do_conj_);
    void perform();
    static void do_all(MatrixDenseView_new<T> * M, char orig_uplo, bool do_conj);
private:
    MatrixDenseView_new<T> * M;
    char orig_uplo = '_';
    bool do_conj = true;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_COMPLETE_DNX_H */
