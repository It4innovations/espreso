
#include "math/wrappers/math.spsolver.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef HAVE_MKL
#ifndef HAVE_PARDISO
#ifndef HAVE_SUITESPARSE

namespace espreso {

struct Solver_External_Representation
{
};

template <typename T, typename I>
const char* DirectSparseSolver<T, I>::name()
{
	return "NONE";
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideFactors()
{
        return false;
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideSC()
{
        return false;
}

template <typename T, typename I>
Solver_Factors DirectSparseSolver<T, I>::factorsSymmetry()
{
	return Solver_Factors::NONE;
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver()
: matrix{}, rows{}, nnzA{}, nnzL{}, memoryL{}, _solver{nullptr}
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(const Matrix_CSR<T> &a)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I> & DirectSparseSolver<T, I>::operator=(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I>::~DirectSparseSolver()
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::commit(const Matrix_CSR<T,I> &a)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::symbolicFactorization(int fixedSuffix)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::numericalFactorization()
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity)
{

}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity)
{

}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixSize()
{
	return 0;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixNnz()
{
	return 0;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getFactorNnz()
{
	return 0;
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T,I> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T,I> &/*U*/, bool /*copyPattern*/, bool /*copyValues*/)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Permutation<I> &/*perm*/)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getSC(Matrix_Dense<T,I> &/*sc*/)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template struct DirectSparseSolver<Matrix_CSR, double, int>;
template struct DirectSparseSolver<Matrix_CSR, std::complex<double>, int>;

}

#endif
#endif
#endif
