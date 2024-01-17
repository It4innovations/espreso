
#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef HAVE_MKL
#ifndef HAVE_PARDISO
#ifndef HAVE_SUITESPARSE

namespace espreso {

template <typename T, typename I, template <typename, typename> class Matrix>
const char* DirectSolver<T, Matrix>::name()
{
	return "NONE";
}

template <typename T, typename I, template <typename, typename> class Matrix>
bool DirectSolver<T, Matrix>::provideFactors()
{
        return false;
}

template <typename T, typename I, template <typename, typename> class Matrix>
bool DirectSolver<T, Matrix>::provideSC()
{
        return false;
}

template <typename T, typename I, template <typename, typename> class Matrix>
DirectSolver<T, Matrix>::DirectSolver()
: matrix{}, rows{}, nnzA{}, nnzL{}, memoryL{}, _solver{nullptr}
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I, template <typename, typename> class Matrix>
DirectSolver<T, Matrix>::DirectSolver(const Matrix<T> &a)
: matrix{}, rows{}, nnzA{}, nnzL{}, memoryL{}, _solver{nullptr}
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, typename I, template <typename, typename> class Matrix>
DirectSolver<T, Matrix>::~DirectSolver()
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::commit(const Matrix<T> &a)
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::commit(SpBLAS<T, Matrix> &spblas)
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::symbolicFactorization(esint fixedSuffix)
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::numericalFactorization()
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution, int sparsity)
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution, int sparsity)
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::getSC(Matrix_Dense<T, I> &sc)
{

}

template <typename T, typename I, template <typename, typename> class Matrix>
void DirectSolver<T, Matrix>::getFactors(Matrix_CSC<T> &L, Matrix_CSC<T> &U, Vector_Dense<int> &p)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template struct DirectSolver<double, Matrix_CSR>;
template struct DirectSolver<std::complex<double>, Matrix_CSR>;

}

#endif
#endif
#endif
