
#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#include <complex>

#ifndef HAVE_MKL
#ifndef HAVE_PARDISO
#ifndef HAVE_SUITESPARSE

namespace espreso {
namespace math {

const char* sparseSolver()
{
	return "NONE";
}

template <>
void initSolver(Matrix_CSR<double> &A)
{

}

template <>
void initSolver(Matrix_CSR<std::complex<double> > &A)
{

}

template <>
void symbolicFactorization(const Matrix_CSR<double> &A, esint fixedSuffix)
{

}

template <>
void symbolicFactorization(const Matrix_CSR<std::complex<double> > &A, esint fixedSuffix)
{

}

template <>
void numericalFactorization(const Matrix_CSR<double> &A)
{

}

template <>
void numericalFactorization(const Matrix_CSR<std::complex<double> > &A)
{

}

template <>
void solve(const Matrix_CSR<double> &A, Vector_Dense<double> &b, Vector_Dense<double> &x, VectorSparsity sparsity)
{

}

template <>
void solve(const Matrix_CSR<double> &A, Matrix_Dense<double> &b, Matrix_Dense<double> &x, VectorSparsity sparsity)
{

}

template <>
void solve(const Matrix_CSR<std::complex<double> > &A, Vector_Dense<std::complex<double> > &b, Vector_Dense<std::complex<double> > &x, VectorSparsity sparsity)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &A, Matrix_Dense<std::complex<double> > &b, Matrix_Dense<std::complex<double> > &x, VectorSparsity sparsity)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <>
void freeSolver(Matrix_CSR<double> &A)
{

}

template <>
void freeSolver(Matrix_CSR<std::complex<double> > &A)
{

}

template <>
void freeFactor(Matrix_CSC<double> &m)
{

}

template <>
void freeFactor(Matrix_CSC<std::complex<double> > &m)
{

}

template <>
SolverInfo getSolverInfo(const Matrix_CSR<double> &m)
{
	return SolverInfo{};
}

template <>
SolverInfo getSolverInfo(const Matrix_CSR<std::complex<double> > &m)
{
	return SolverInfo{};
}

bool provideFactors()
{
	return false;
}

bool provideSC()
{
	return false;
}

template <>
void computeSC(const Matrix_CSR<double> &m, Matrix_Dense<double> &sc)
{

}

template <>
void computeSC(const Matrix_CSR<std::complex<double> > &m, Matrix_Dense<std::complex<double> > &sc)
{

}

template <>
void getFactors(const Matrix_CSR<double> &m, Matrix_CSC<double> &L, Matrix_CSC<double> &U, Vector_Dense<int> &p)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <>
void getFactors(const Matrix_CSR<std::complex<double> > &m, Matrix_CSC<std::complex<double> > &L, Matrix_CSC<std::complex<double> > &U, Vector_Dense<int> &p)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}


}
}

#endif
#endif
#endif
