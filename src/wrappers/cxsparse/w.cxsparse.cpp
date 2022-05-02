
#include "w.cxsparse.h"

#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_CXSPARSE
#ifdef USE_SOLVER_CXSPARSE

namespace espreso {

struct Matrix_CSR_Solver { };

namespace math {


template <>
void initSolver(Matrix_CSR<double> &x)
{
	x._solver = new Matrix_CSR_Solver();
}

template <>
void initSolver(Matrix_CSR<std::complex<double> > &x)
{
	x._solver = new Matrix_CSR_Solver();
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

}

template <>
void solve(const Matrix_CSR<std::complex<double> > &x, Matrix_Dense<std::complex<double> > &rhs, Matrix_Dense<std::complex<double> > &solution)
{

}

template <>
void freeSolver(Matrix_CSR<double> &x)
{
	delete x._solver;
}

template <>
void freeSolver(Matrix_CSR<std::complex<double> > &x)
{
	delete x._solver;
}

}
}

#endif
#endif

