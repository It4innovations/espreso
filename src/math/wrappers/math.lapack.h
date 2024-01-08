
#ifndef SRC_MATH_WRAPPERS_MATH_LAPACK_H_
#define SRC_MATH_WRAPPERS_MATH_LAPACK_H_

#include "math/primitives/matrix_dense.h"

namespace espreso {
namespace math {
namespace lapack {

template <typename T>
void solve(Matrix_Dense<T> &A, Matrix_Dense<T> &rhs);

}
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_LAPACK_H_ */
