
#include "matrix.csr.feti.h"
#include "matrix.dense.feti.h"

using namespace espreso;

MatrixCSRFETI::MatrixCSRFETI()
{

}

MatrixCSRFETI::MatrixCSRFETI(const MatrixCSRFETI &other)
{
	deepCopy(&other);
}

MatrixCSRFETI& MatrixCSRFETI::operator=(const MatrixCSRFETI &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixCSRFETI* MatrixCSRFETI::copy()
{
	return new MatrixCSRFETI();
}

void MatrixCSRFETI::multiply(MatrixCSRFETI &A, MatrixCSRFETI &B, bool transposeA)
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->multiply(A[d], B[d], transposeA);
	}
}

void MatrixCSRFETI::solve(const MatrixDenseFETI &rhs, MatrixDenseFETI &solution)
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->solve(rhs[d], solution[d]);
	}
}

void MatrixCSRFETI::removeLower(MatrixType type)
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->removeLower(type);
	}
}

Matrix* MatrixCSRFETI::create()
{
	return new MatrixCSR();
}





