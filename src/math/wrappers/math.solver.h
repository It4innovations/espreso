
#ifndef SRC_MATH_WRAPPERS_MATH_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SOLVER_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

namespace espreso {
namespace math {

	template <typename T> void initSolver(Matrix_CSR<T> &x);

	template <typename T> void symbolicFactorization(const Matrix_CSR<T> &x);
	template <typename T> void numericalFactorization(const Matrix_CSR<T> &x);
	template <typename T> void solve(const Matrix_CSR<T> &x, Vector_Dense<T> &rhs, Vector_Dense<T> &solution);
	template <typename T> void solve(const Matrix_CSR<T> &x, Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution);

	template <typename T> void freeSolver(Matrix_CSR<T> &x);
}
}



#endif /* SRC_MATH_WRAPPERS_MATH_SOLVER_H_ */
