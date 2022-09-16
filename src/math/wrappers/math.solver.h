
#ifndef SRC_MATH_WRAPPERS_MATH_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SOLVER_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"

#include <cstddef>

namespace espreso {
namespace math {

	enum VectorSparsity {
		DENSE           = 0,
		SPARSE_RHS      = 1 << 0,
		SPARSE_SOLUTION = 1 << 1
	};

	struct SolverInfo {
		size_t rows, nnzA, nnzL;
		size_t memoryL;
	};

	const char* sparseSolver();

	template <typename T> void initSolver(Matrix_CSR<T> &m);

	template <typename T> void symbolicFactorization(const Matrix_CSR<T> &m, esint fixedSuffix = 0); // do not permute suffix
	template <typename T> void numericalFactorization(const Matrix_CSR<T> &m);
	template <typename T> void solve(const Matrix_CSR<T> &m, Vector_Dense<T> &rhs, Vector_Dense<T> &solution, VectorSparsity sparsity = VectorSparsity::DENSE);
	template <typename T> void solve(const Matrix_CSR<T> &m, Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution, VectorSparsity sparsity = VectorSparsity::DENSE);

	bool provideFactors();
	bool provideSC();

	template <typename T> SolverInfo getSolverInfo(const Matrix_CSR<T> &m);
	template <typename T> void getFactors(const Matrix_CSR<T> &m, Matrix_CSC<T> &L, Matrix_CSC<T> &U, Vector_Dense<int> &p);
	template <typename T> void computeSC(const Matrix_CSR<T> &m, Matrix_Dense<T> &sc);

	template <typename T> void freeSolver(Matrix_CSR<T> &m);
	template <typename T> void freeFactor(Matrix_CSC<T> &m);

inline VectorSparsity operator|(const VectorSparsity &s1, const VectorSparsity &s2) { return (VectorSparsity)((int)s1 | (int)s2); }
}
}

#endif /* SRC_MATH_WRAPPERS_MATH_SOLVER_H_ */

