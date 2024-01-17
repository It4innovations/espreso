
#ifndef SRC_MATH_WRAPPERS_MATH_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SOLVER_H_

#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"
#include "math.spblas.h"

#include <cstddef>

namespace espreso {

struct Matrix_Solver_External_Representation;

template <template <typename, typename> class Matrix, typename T, typename I = int>
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
	DirectSolver(const Matrix<T, I> &a);
	~DirectSolver();

	void commit(const Matrix<T, I> &a);
	void commit(SpBLAS<Matrix, T, I> &spblas);

	void symbolicFactorization(int fixedSuffix = 0); // do not permute suffix
	void numericalFactorization();

	void solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity = VectorSparsity::DENSE);
	void solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity = VectorSparsity::DENSE);

	void getFactors(Matrix_CSC<T> &L, Matrix_CSC<T> &U, Vector_Dense<int> &p);
	void getSC(Matrix_Dense<T, I> &sc);

	const Matrix<T, I> *matrix;

	size_t rows, nnzA, nnzL;
	size_t memoryL;

private:
	Matrix_Solver_External_Representation *_solver;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_SOLVER_H_ */

