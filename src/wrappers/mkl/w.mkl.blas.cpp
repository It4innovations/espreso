
#include "math/math.h"
#include "math/math.h"

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
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_Dense<double> &a, const double &beta, const Vector_Dense<double> &x)
{
#ifdef HAVE_MKL
	cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, a.nrows, a.ncols, alpha, a.vals, a.ncols, x.vals, 1, beta, y.vals, 1);
#endif
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_Dense<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	cblas_zgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, a.nrows, a.ncols, &alpha, a.vals, a.ncols, x.vals, 1, &beta, y.vals, 1);
#endif
}

}
}


