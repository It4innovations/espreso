
#include "math/math.h"
#include "math/matrix.h"
#include "esinfo/eslog.h"

using namespace espreso;

#ifndef HAVE_PARDISO
#ifdef HAVE_MKL
#include "w.mkl.h"
#include "mkl_pardiso.h"

static void call(esint phase, MatrixType type, MATH::CSRHandler *A, esint nrhs, double *rhs, double *solution)
{
	esint nrows, ncols, nnz, *rows, *cols;
	double *vals;
	A->inner->info(nrows, ncols, nnz, rows, cols, vals);

	switch (type) {
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		A->inner->mtype = 2; break;
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
		A->inner->mtype = -2; break;
	case MatrixType::REAL_UNSYMMETRIC:
		A->inner->mtype = 1; break;
	}

	A->inner->phase = phase;
	pardiso(A->inner->pt,
			&A->inner->maxfct, &A->inner->mnum,
			&A->inner->mtype, &A->inner->phase,
			&nrows, vals, rows, cols,
			A->inner->perm,
			&nrhs,
			A->inner->iparm, &A->inner->msglvl,
			rhs, solution,
			&A->inner->error);

	switch (A->inner->error) {
	case   0: break;
	case  -1: eslog::error("MKL PARDISO: input inconsistent.\n"); break;
	case  -2: eslog::error("MKL PARDISO: not enough memory.\n"); break;
	case  -3: eslog::error("MKL PARDISO: reordering problem.\n"); break;
	case  -4: eslog::error("MKL PARDISO: zero pivot, numerical factorization or iterative refinement problem.\n"); break;
	case  -5: eslog::error("MKL PARDISO: unclassified (internal) error.\n"); break;
	case  -6: eslog::error("MKL PARDISO: reordering failed.\n"); break;
	case  -7: eslog::error("MKL PARDISO: diagonal matrix is singular.\n"); break;
	case  -8: eslog::error("MKL PARDISO: 32-bit integer overflow problem.\n"); break;
	case  -9: eslog::error("MKL PARDISO: not enough memory for OOC.\n"); break;
	case -10: eslog::error("MKL PARDISO: error opening OOC files.\n"); break;
	case -11: eslog::error("MKL PARDISO: read/write error with OOC files.\n"); break;
	}
}
#endif

void MATH::CSRMatFactorizeSymbolic(MatrixType type, CSRHandler *A)
{
#ifdef HAVE_MKL
	call(11, type, A, 1, NULL, NULL);
	A->inner->called = 1;
#endif
}

void MATH::CSRMatFactorizeNumeric(MatrixType type, CSRHandler *A)
{
#ifdef HAVE_MKL
	if (A->inner->called == 0) {
		call(12, type, A, 1, NULL, NULL);
	}
	if (A->inner->called == 1) {
		call(22, type, A, 1, NULL, NULL);
	}
	A->inner->called = 2;
#endif
}

void MATH::CSRMatFactorize(MatrixType type, CSRHandler *A)
{
#ifdef HAVE_MKL
	call(12, type, A, 1, NULL, NULL);
	A->inner->called = 2;
#endif
}

void MATH::CSRMatSolve(MatrixType type, CSRHandler *A, esint nrhs, double *rhs, double *solution)
{
#ifdef HAVE_MKL
	if (A->inner->called < 2) {
		CSRMatFactorize(type, A);
	}
	call(33, type, A, nrhs, rhs, solution);
#endif
}

void MATH::CSRMatClearFactors(CSRHandler *A)
{
#ifdef HAVE_MKL
	if (A->inner->called) {
		call(-1, MatrixType::REAL_UNSYMMETRIC, A, 1, NULL, NULL);
	}
#endif
}

#endif // HAVE_PARDISO
