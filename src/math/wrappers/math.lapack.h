
#ifndef SRC_MATH_WRAPPERS_MATH_LAPACK_H_
#define SRC_MATH_WRAPPERS_MATH_LAPACK_H_

#include "math/primitives/matrix_dense.h"
#include "math/primitives/vector_dense.h"

namespace espreso {
namespace math {
namespace lapack {

template <typename T, typename I>
void solve_sym_upper(Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &rhs);

template <typename T, typename I>
void solve_general(Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &rhs);

}
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_LAPACK_H_ */
