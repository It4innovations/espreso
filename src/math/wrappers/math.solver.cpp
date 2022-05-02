
#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#ifndef HAVE_MKL
#ifndef HAVE_PARDISO

namespace espreso {
namespace math {

template <>
void initSolver(const Matrix_CSR<double> &x)
{

}

template <>
void initSolver(const Matrix_CSR<std::complex<double> > &x)
{

}

template <>
void symbolicFactorization(const Matrix_CSR<double> &x)
{

}

template <>
void symbolicFactorization(const Matrix_CSR<std::complex<double> > &x)
{

}

template <>
void numericalFactorization(const Matrix_CSR<double> &x)
{

}

template <>
void numericalFactorization(const Matrix_CSR<std::complex<double> > &x)
{

}

template <>
void solve(const Matrix_CSR<double> &x, Vector_Dense<double> &rhs, Vector_Dense<double> &solution)
{

}

template <>
void solve(const Matrix_CSR<double> &x, Matrix_Dense<double> &rhs, Matrix_Dense<double> &solution)
{

}

template <>
void solve(const Matrix_CSR<std::complex<double> > &x, Vector_Dense<std::complex<double> > &rhs, Vector_Dense<std::complex<double> > &solution)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &x, Matrix_Dense<std::complex<double> > &rhs, Matrix_Dense<std::complex<double> > &solution)
{
	eslog::error("calling of empty sparse solver wrapper.\n");
}

template <>
void freeSolver(const Matrix_CSR<double> &x)
{

}

template <>
void freeSolver(const Matrix_CSR<std::complex<double> > &x)
{

}

}
}

#endif
#endif
