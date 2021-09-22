
#include "mklpdsssolver.h"
#include "mklpdsssystem.h"

void espreso::setDirichlet(Matrix_Distributed<Matrix_CSR, double> &A, Vector_Distributed<Vector_Dense, double> &b, const Vector_Sparse<double> &dirichlet, const DOFsDistribution &distribution)
{
	_setDirichlet<double>(A, b, dirichlet, distribution);
}

void espreso::setDirichlet(Matrix_Distributed<Matrix_CSR, std::complex<double>> &A, Vector_Distributed<Vector_Dense, std::complex<double>> &b, const Vector_Sparse<std::complex<double>> &dirichlet, const DOFsDistribution &distribution)
{
	_setDirichlet<std::complex<double>>(A, b, dirichlet, distribution);
}

