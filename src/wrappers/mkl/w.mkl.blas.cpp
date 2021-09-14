
#include "math/math.h"
#include "math2/math2.h"

#ifdef HAVE_MKL
#include "mkl_blas.h"
#include "mkl_cblas.h"
#endif

using namespace espreso;

namespace espreso {
namespace math {

template <>
void copy(const esint size, double *x, const esint incX, const double *y, const esint incY)
{
#ifdef HAVE_MKL
	cblas_dcopy(size, y, incY, x, incX);
#endif
}

template <>
void copy(const esint size, std::complex<double> *x, const esint incX, const std::complex<double> *y, const esint incY)
{
#ifdef HAVE_MKL
	cblas_zcopy(size, y, incY, x, incX);
#endif
}

template <>
void scale(const esint size, const double &alpha, double *x, const esint incX)
{
#ifdef HAVE_MKL
	cblas_dscal(size, alpha, x, incX);
#endif
}

template <>
void scale(const esint size, const std::complex<double> &alpha, std::complex<double> *x, const esint incX)
{
#ifdef HAVE_MKL
	cblas_zscal(size, &alpha, x, incX);
#endif
}

template <>
void add(const esint size, double *x, const esint incX, const double &alpha, const double *y, const esint incY)
{
#ifdef HAVE_MKL
	cblas_daxpy(size, alpha, y, incY, x, incX);
#endif
}

template <>
void add(const esint size, std::complex<double> *x, const esint incX, const std::complex<double> &alpha, const std::complex<double> *y, const esint incY)
{
#ifdef HAVE_MKL
	cblas_zaxpy(size, &alpha, y, incY, x, incX);
#endif
}

template <>
double dot(const esint size, const double *x, const esint incX, const double *y, const esint incY)
{
	double dot = 0;
#ifdef HAVE_MKL
	dot = cblas_ddot(size, x, incX, y, incY);
#endif
	return dot;
}

template <>
std::complex<double> dot(const esint size, const std::complex<double> *x, const esint incX, const std::complex<double> *y, const esint incY)
{
	std::complex<double> dot = 0;
#ifdef HAVE_MKL
	cblas_cdotu_sub(size, x, incX, y, incY, &dot);
#endif
	return dot;
}

template <>
double norm(const esint size, const double *x, const esint incX)
{
	double norm = 0;
#ifdef HAVE_MKL
	norm = dnrm2(&size, x, &incX);
#endif
	return norm;
}

template <>
void apply(Vector_Dense<double> &y, const double &alpha, Matrix_Dense<double> &a, const double &beta, const Vector_Dense<double> &x)
{
#ifdef HAVE_MKL
	cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, a.nrows, a.ncols, alpha, a.vals, a.ncols, x.vals, 1, beta, y.vals, 1);
#endif
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, Matrix_Dense<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	cblas_zgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, a.nrows, a.ncols, &alpha, a.vals, a.ncols, x.vals, 1, &beta, y.vals, 1);
#endif
}

}
}


/// OLD ///

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

