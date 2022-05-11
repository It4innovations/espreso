
#ifndef SRC_MATH_WRAPPERS_MATH_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SOLVER_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

namespace espreso {
namespace math {

	template <typename T> void initSolver(Matrix_CSR<T> &m);
	template <typename T> void restrictToSurface(Matrix_CSR<T> &m, esint surfaceSize);

	template <typename T> void symbolicFactorization(const Matrix_CSR<T> &m);
	template <typename T> void numericalFactorization(const Matrix_CSR<T> &m);
	template <typename T> void solve(const Matrix_CSR<T> &m, Vector_Dense<T> &rhs, Vector_Dense<T> &solution);
	template <typename T> void solve(const Matrix_CSR<T> &m, Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution);

	template <typename T> void computeSC(const Matrix_CSR<T> &m, Matrix_Dense<T> &sc);

	template <typename T> void freeSolver(Matrix_CSR<T> &m);
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_SOLVER_H_ */
