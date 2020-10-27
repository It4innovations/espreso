
#include "math/math.h"
#include "esinfo/eslog.h"

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

void MATH::DenseMinGeneralizedEigenVectors(esint msize, double *A, double *B, esint n, double *lambdas, double *vectors)
{
#ifdef HAVE_MKL
	esint nn;
	esint *fail = new esint[msize];
	double *tmpvec = new double[msize * msize];
	LAPACKE_dsygvx(LAPACK_ROW_MAJOR, 1, 'V', 'I', 'U', msize, A, msize, B, msize, 0, 0, 1, n, 0, &nn, lambdas, tmpvec, msize, fail);
	for (esint r = 0; r < msize; ++r) {
		for (esint c = 0; c < n; ++c) {
			vectors[r * n + c] = tmpvec[r * msize + c];
		}
	}
	delete[] fail;
	delete[] tmpvec;
#endif
}

void MATH::upDense3x3EigenValuesEigenVectors(double *A, double *W, double *Z)
{
#ifdef HAVE_MKL
	for(int i = 0; i < 9; ++i){
		Z[i] = A[i];
	}
	esint info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', 3, Z, 3, W );
	
	if(info > 0){
		eslog::globalerror("ESPRESO internal error: upDense3x3EigenValuesEigenVectors failed to compued eigen-values or eigen-vectors.\n");
	}
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
