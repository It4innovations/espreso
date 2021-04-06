
#include "math/math.h"

#ifdef HAVE_MKL
#include "mkl_lapacke.h"
#endif

using namespace espreso;

void MATH::upDense3x3EigenValues(double *mVals, double *eigenValues)
{
#ifdef HAVE_MKL
	double T[6] = { mVals[0], mVals[1], mVals[2], mVals[3], mVals[4], mVals[5] };
	LAPACKE_dsterf(3, T, T + 3);
	eigenValues[0] = T[2];
	eigenValues[1] = T[1];
	eigenValues[2] = T[0];
#endif
}

void MATH::DenseMatDenseMatRowMajorSystemSolve(int nra, int ncb, double *a, double *b)
{
	// B = A\B
#ifdef HAVE_MKL
	//lapack_int LAPACKE_dgesv (int matrix_layout, lapack_int n, lapack_int nrhs, double * a, lapack_int lda, lapack_int * ipiv, double * b, lapack_int ldb)
	esint ipiv[nra];
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, nra, ncb, a, nra, ipiv, b, ncb);
#endif
}
