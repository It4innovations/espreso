
#include "math/math.h"

#ifdef HAVE_MKL
#include "mkl_blas.h"
#include "mkl_cblas.h"
#endif

using namespace espreso;

void MATH::vecScale(esint size, float alpha, float *vVals)
{
#ifdef HAVE_MKL
	esint incr = 1;
	cblas_sscal(size, alpha, vVals, incr);
#endif
}

void MATH::vecScale(esint size, double alpha, double *vVals)
{
#ifdef HAVE_MKL
	esint incr = 1;
	cblas_dscal(size, alpha, vVals, incr);
#endif
}

double MATH::vecDot(esint size, double *vVals)
{
	double dot = 0;
#ifdef HAVE_MKL
	esint incr = 1;
	dot = cblas_ddot(size, vVals, incr, vVals, incr);
#endif
	return dot;
}

double MATH::vecDot(esint size, double *a, double *b)
{
	double dot = 0;
#ifdef HAVE_MKL
	esint incr = 1;
	dot = cblas_ddot(size, a, incr, b, incr);
#endif
	return dot;
}

void MATH::vecAdd(esint size, double *result, double alpha, double *other)
{
#ifdef HAVE_MKL
	esint incr = 1;
	cblas_daxpy(size, alpha, other, incr, result, incr);
#endif
}

void MATH::vecAddSparse(esint size, double *result, double alpha, esint *indices, double *other)
{
#ifdef HAVE_MKL
	cblas_daxpyi(size, alpha, other, indices, result);
#endif
}

void MATH::vecDenseToSparse(esint size, esint *indices, double *sparse, double *dense)
{
#ifdef HAVE_MKL
	cblas_dgthr(size, dense, sparse, indices);
#endif
}

double MATH::vecNorm(esint size, float *vVals)
{
	double norm = 0;
#ifdef HAVE_MKL
	esint incr = 1;
	norm = snrm2(&size, vVals, &incr);
#endif
	return norm;
}

double MATH::vecNorm(esint size, double *vVals)
{
	double norm = 0;
#ifdef HAVE_MKL
	esint incr = 1;
	norm = dnrm2(&size, vVals, &incr);
#endif
	return norm;
}

esint MATH::vecNormMaxIndex(esint size, float *vVals)
{
	esint index = -1;
#ifdef HAVE_MKL
	esint incr = 1;
	index = cblas_isamax(size, vVals, incr);
#endif
	return index;
}

esint MATH::vecNormMaxIndex(esint size, double *vVals)
{
	esint index = -1;
#ifdef HAVE_MKL
	esint incr = 1;
	index = cblas_idamax(size, vVals, incr);
#endif
	return index;
}

void MATH::DenseMatDenseMatRowMajorProduct(
			double alpha, bool transposeA, esint aRows, esint aCols, double* aVals,
			bool transposeB, esint bRows, esint bCols, double* bVals,
			double beta, double* cVals)
{
#ifdef HAVE_MKL
	cblas_dgemm(
		CblasRowMajor,
		transposeA ? CblasTrans : CblasNoTrans,
		transposeB ? CblasTrans : CblasNoTrans,
		transposeA ? aCols : aRows,
		transposeB ? bRows : bCols,
		transposeA ? aRows : aCols,
		alpha,
		aVals, aCols,
		bVals, bCols,
		beta,
		cVals, transposeB ? bRows : bCols);
#endif
}

