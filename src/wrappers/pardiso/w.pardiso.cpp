
#include "math/math.h"
#include "math/math.h"
#include "math/matrix.h"
#include "esinfo/eslog.h"
#include "esinfo/envinfo.h"

#ifdef HAVE_MKL
#include "wrappers/mkl/w.mkl.h"
#endif

using namespace espreso;

#ifdef HAVE_PARDISO

extern "C" {
void pardisoinit(void*, int*, int*, int*, double*, int*);
void pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);
void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
void pardiso_chkvec(int*, int*, double*, int*);
void pardiso_printstats(int*, int*, double*, int*, int*, int*, double*, int*);
}

static inline void setMatrixType(MatrixType type, MATH::CSRHandler *A)
{
	switch (type) {
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		A->inner->mtype = 2; break;
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
		A->inner->mtype = -2; break;
	case MatrixType::REAL_UNSYMMETRIC:
		A->inner->mtype = 1; break;
	}
}

static void call(esint phase, MatrixType type, MATH::CSRHandler *A, esint nrhs, double *rhs, double *solution)
{
	esint nrows, ncols, nnz, *rows, *cols;
	double *vals;
	A->inner->info(nrows, ncols, nnz, rows, cols, vals);

	setMatrixType(type, A);

	A->inner->phase = phase;
	pardiso(A->inner->pt,
		&A->inner->maxfct, &A->inner->mnum,
		&A->inner->mtype, &A->inner->phase,
		&nrows, vals, rows, cols,
		A->inner->perm,
		&nrhs,
		A->inner->iparm, &A->inner->msglvl,
		rhs, solution,
		&A->inner->error,
		A->inner->dparm);

	switch (A->inner->error) {
	case   0: break;
	case  -1: eslog::error("PARDISO: input inconsistent.\n"); break;
	case  -2: eslog::error("PARDISO: not enough memory.\n"); break;
	case  -3: eslog::error("PARDISO: reordering problem.\n"); break;
	case  -4: eslog::error("PARDISO: zero pivot, numerical factorization or iterative refinement problem.\n"); break;
	case  -5: eslog::error("PARDISO: unclassified (internal) error.\n"); break;
	case  -6: eslog::error("PARDISO: reordering failed.\n"); break;
	case  -7: eslog::error("PARDISO: diagonal matrix is singular.\n"); break;
	case  -8: eslog::error("PARDISO: 32-bit integer overflow problem.\n"); break;
	case  -9: eslog::error("PARDISO: not enough memory for OOC.\n"); break;
	case -10: eslog::error("PARDISO: error opening OOC files.\n"); break;
	case -11: eslog::error("PARDISO: read/write error with OOC files.\n"); break;
	}
}

static void callInit(MatrixType type, MATH::CSRHandler *A)
{
	setMatrixType(type, A);

	int solver = 0;   // direct solver
	pardisoinit(A->inner->pt,
		&A->inner->mtype,
		&solver,
		A->inner->iparm,
		A->inner->dparm,
		&A->inner->error);

	A->inner->iparm[3 - 1] = info::env::OMP_NUM_THREADS;

	switch (A->inner->error) {
	case   0: break;
	case -10: eslog::error("PARDISO: No license file pardiso.lic found.\n"); break;
	case -11: eslog::error("PARDISO: License is expired.\n"); break;
	case -12: eslog::error("PARDISO: Wrong username or hostname.\n"); break;
	}
}

void MATH::CSRMatFactorizeSymbolic(MatrixType type, CSRHandler *A)
{
	if (A->inner->called < 1) {
		callInit(type, A);
	}

	call(11, type, A, 1, NULL, NULL);
	A->inner->called = 1;
}

void MATH::CSRMatFactorizeNumeric(MatrixType type, CSRHandler *A)
{
	if (A->inner->called < 1) {
		CSRMatFactorizeSymbolic(type, A);
	}

	call(22, type, A, 1, NULL, NULL);
	A->inner->called = 2;
}

void MATH::CSRMatFactorize(MatrixType type, CSRHandler *A)
{
	CSRMatFactorizeNumeric(type, A);
}

void MATH::CSRMatSolve(MatrixType type, CSRHandler *A, esint nrhs, double *rhs, double *solution)
{
	if (A->inner->called < 2) {
		CSRMatFactorizeNumeric(type, A);
	}
	call(33, type, A, nrhs, rhs, solution);
}

void MATH::CSRMatClearFactors(CSRHandler *A)
{
	if (A->inner->called) {
		call(-1, MatrixType::REAL_UNSYMMETRIC, A, 1, NULL, NULL);
		A->inner->called = 0;
	}
}

#endif
