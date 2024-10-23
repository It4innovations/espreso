
#ifndef SRC_MATH_WRAPPERS_MATH_SPSOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SPSOLVER_H_

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

enum struct Solver_Factors: int {
	NONE,
	NONSYMMETRIC_BOTH,
	HERMITIAN_UPPER, // for non-complex matrices, hermitian and symmetric is equivalent
	HERMITIAN_LOWER
};

template<typename T, typename I>
struct Solver_External_Representation;

template <typename T, typename I = int>
struct DirectSparseSolver {
//    struct VectorSparsity {
//        static const int DENSE           = 0;
//        static const int SPARSE_RHS      = 1 << 0;
//        static const int SPARSE_SOLUTION = 2 << 0;
//    };

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

    void symbolicFactorization();
    void numericalFactorization();

    void solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution);
    void solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution);

    void solveForward (Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution);
    void solveDiagonal(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution);
    void solveBackward(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution);
    void solveForward (Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution);
    void solveDiagonal(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution);
    void solveBackward(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution);


    I getMatrixSize();
    I getMatrixNnz();
    I getFactorNnz();

    template<typename A>
    void getFactorL(Matrix_CSR<T,I,A> &L, bool copyPattern = true, bool copyValues = true);
    template<typename A>
    void getFactorU(Matrix_CSR<T,I,A> &U, bool copyPattern = true, bool copyValues = true);
    void getPermutation(Permutation<I> &perm);
    void getPermutation(Vector_Dense<I> &perm);
    void getSC(Matrix_Dense<T,I> &sc);

private:
    std::unique_ptr<Solver_External_Representation<T,I>> ext;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_SPSOLVER_H_ */

