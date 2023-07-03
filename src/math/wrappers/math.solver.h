
#ifndef SRC_MATH_WRAPPERS_MATH_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SOLVER_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"

#include <cstddef>

namespace espreso {

struct Matrix_Solver_External_Representation;

template <typename T, template <typename> class Matrix>
struct DirectSolver {
	struct VectorSparsity {
		static const int DENSE           = 0;
		static const int SPARSE_RHS      = 1 << 0;
		static const int SPARSE_SOLUTION = 2 << 0;
	};

	static const char* name();
	static bool provideFactors();
	static bool provideSC();

	DirectSolver();
	DirectSolver(const Matrix<T> &a);
	~DirectSolver();

	void commit(const Matrix<T> &a);

	void symbolicFactorization(int fixedSuffix = 0); // do not permute suffix
	void numericalFactorization();

	void solve(Vector_Dense<T> &rhs, Vector_Dense<T> &solution, int sparsity = VectorSparsity::DENSE);
	void solve(Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution, int sparsity = VectorSparsity::DENSE);

	void getFactors(Matrix_CSC<T> &L, Matrix_CSC<T> &U, Vector_Dense<int> &p);
	void getSC(Matrix_Dense<T> &sc);

	const Matrix<T> *matrix;

	size_t rows, nnzA, nnzL;
	size_t memoryL;

private:
	Matrix_Solver_External_Representation *_solver;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_SOLVER_H_ */

