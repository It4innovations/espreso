
#ifndef SRC_MATH_WRAPPERS_MATH_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SOLVER_H_

#include "basis/containers/allocators.h"
#include "math/primitives/matrix_info.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_csc.h"
#include "math/primitives/permutation.h"
#include "math.spblas.h"

#include <cstddef>

namespace espreso {

template<typename T, typename I>
struct Solver_External_Representation;

template <typename T, typename I = int>
struct DirectSparseSolver {
	struct VectorSparsity {
		static const int DENSE           = 0;
		static const int SPARSE_RHS      = 1 << 0;
		static const int SPARSE_SOLUTION = 2 << 0;
	};

	static const char* name();
	static bool provideFactors();
	static bool provideSC();
	static Solver_Factors factorsSymmetry();

	DirectSparseSolver();
	DirectSparseSolver(const Matrix_CSR<T, I> &a);
	DirectSparseSolver(const DirectSparseSolver &other) = delete;
	DirectSparseSolver(DirectSparseSolver &&other);
	DirectSparseSolver & operator=(const DirectSparseSolver &other) = delete;
	DirectSparseSolver & operator=(DirectSparseSolver &&other);
	~DirectSparseSolver();

	void commit(const Matrix_CSR<T, I> &a);

	void symbolicFactorization(int fixedSuffix = 0); // do not permute suffix
	void numericalFactorization();

	void solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity = VectorSparsity::DENSE);
	void solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity = VectorSparsity::DENSE);

	I getMatrixSize();
	I getMatrixNnz();
	I getFactorNnz();

	void getFactorL(Matrix_CSR<T,I> &L, bool copyPattern = true, bool copyValues = true);
	void getFactorU(Matrix_CSR<T,I> &U, bool copyPattern = true, bool copyValues = true);
	void getPermutation(Permutation<I> &perm);
	void getSC(Matrix_Dense<T,I> &sc);

private:
	std::unique_ptr<Solver_External_Representation<T,I>> ext;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_SOLVER_H_ */

