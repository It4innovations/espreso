
#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef MKL
#ifndef PARDISO
#ifndef SUITESPARSE

namespace espreso {

template <typename T, template <typename> class Matrix>
const char* DirectSolver<T, Matrix>::name()
{
	return "NONE";
}

template <typename T, template <typename> class Matrix>
bool DirectSolver<T, Matrix>::provideFactors()
{
        return false;
}

template <typename T, template <typename> class Matrix>
bool DirectSolver<T, Matrix>::provideSC()
{
        return false;
}

template <typename T, template <typename> class Matrix>
DirectSolver<T, Matrix>::DirectSolver()
: matrix{}, rows{}, nnzA{}, nnzL{}, memoryL{}, _solver{nullptr}
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, template <typename> class Matrix>
DirectSolver<T, Matrix>::DirectSolver(const Matrix<T> &a)
: matrix{}, rows{}, nnzA{}, nnzL{}, memoryL{}, _solver{nullptr}
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <typename T, template <typename> class Matrix>
DirectSolver<T, Matrix>::~DirectSolver()
{

}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::commit(const Matrix<T> &a)
{

}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::commit(SpBLAS<T, Matrix> &spblas)
{

}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::symbolicFactorization(esint fixedSuffix)
{

}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::numericalFactorization()
{

}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Vector_Dense<T> &rhs, Vector_Dense<T> &solution, int sparsity)
{

}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::solve(Matrix_Dense<T> &rhs, Matrix_Dense<T> &solution, int sparsity)
{

}

template <typename T, template <typename> class Matrix>
void DirectSolver<T, Matrix>::getSC(Matrix_Dense<T> &sc)
{

}

template <typename T, template <typename> class Matrix>
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
