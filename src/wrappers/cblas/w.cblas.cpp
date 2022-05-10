
#include "esinfo/eslog.h"
#include "math/wrappers/math.blas.h"

#include <complex>

#ifdef HAVE_CBLAS
#ifdef USE_BLAS_CBLAS
#include "cblas.h"

namespace espreso {
namespace math {

template <>
void copy(const esint size, double *x, const esint incX, const double *y, const esint incY)
{
	cblas_dcopy(size, y, incY, x, incX);
}

template <>
void copy(const esint size, std::complex<double> *x, const esint incX, const std::complex<double> *y, const esint incY)
{
	cblas_zcopy(size, y, incY, x, incX);
}

template <>
void scale(const esint size, const double &alpha, double *x, const esint incX)
{
	cblas_dscal(size, alpha, x, incX);
}

template <>
void scale(const esint size, const std::complex<double> &alpha, std::complex<double> *x, const esint incX)
{
	cblas_zscal(size, &alpha, x, incX);
}

template <>
void add(const esint size, double *x, const esint incX, const double &alpha, const double *y, const esint incY)
{
	cblas_daxpy(size, alpha, y, incY, x, incX);
}

template <>
void add(const esint size, std::complex<double> *x, const esint incX, const std::complex<double> &alpha, const std::complex<double> *y, const esint incY)
{
	cblas_zaxpy(size, &alpha, y, incY, x, incX);
}

template <>
double dot(const esint size, const double *x, const esint incX, const double *y, const esint incY)
{
	return cblas_ddot(size, x, incX, y, incY);
}

template <>
std::complex<double> dot(const esint size, const std::complex<double> *x, const esint incX, const std::complex<double> *y, const esint incY)
{
	std::complex<double> dot = 0;
	eslog::error("not implemented BLAS routine.\n");
//	return cblas_cdotu_sub(size, x, incX, y, incY, &dot);
	return dot;
}

template <>
double norm(const esint size, const double *x, const esint incX)
{
	return cblas_dnrm2(size, x, incX);
}

template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_Dense<double> &a, const double &beta, const Vector_Dense<double> &x)
{
	if (a.submatrix[0].step != 1 || a.submatrix[1].step != 1) {
		eslog::error("slice is incompatible with apply.");
	}
	Slice rows = a.submatrix[0], cols = a.submatrix[1];
	rows.evaluate(a.nrows); cols.evaluate(a.ncols);
	if (a.shape == Matrix_Shape::FULL) {
		double *vals = a.vals + rows.start * a.ncols + cols.start;
		cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, rows.end - rows.start, cols.end - cols.start, alpha, vals, a.ncols, x.vals, 1, beta, y.vals, 1);
	} else {
		if (rows.start != cols.start || rows.end != cols.end) {
			eslog::error("cannot apply non-square slice to a packed matrix.\n");
		}
		double *vals = a.vals + rows.start * a.ncols + cols.start - (rows.start) * (rows.end - rows.start + 1) / 2;
		CBLAS_UPLO uplo = a.shape == Matrix_Shape::UPPER ? CBLAS_UPLO::CblasUpper : CBLAS_UPLO::CblasLower;
		cblas_dspmv(CBLAS_LAYOUT::CblasRowMajor, uplo, rows.end - rows.start, alpha, vals, x.vals, 1, beta, y.vals, 1);
	}
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_Dense<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	if (a.submatrix[0].step != 1 || a.submatrix[1].step != 1) {
		eslog::error("slice is incompatible with apply.");
	}
	Slice rows = a.submatrix[0], cols = a.submatrix[1];
	rows.evaluate(a.nrows); cols.evaluate(a.ncols);
	if (a.shape == Matrix_Shape::FULL) {
		std::complex<double> *vals = a.vals + rows.start * a.ncols + cols.start;
		cblas_zgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, a.nrows, a.ncols, &alpha, vals, a.ncols, x.vals, 1, &beta, y.vals, 1);
	} else {
		eslog::error("not implemented BLAS routine.\n");
	}
}

template <>
void applyT(Vector_Dense<double> &y, const double &alpha, const Matrix_Dense<double> &a, const double &beta, const Vector_Dense<double> &x)
{
	if (a.submatrix[0].step != 1 || a.submatrix[1].step != 1) {
		eslog::error("slice is incompatible with apply.");
	}
	Slice rows = a.submatrix[0], cols = a.submatrix[1];
	rows.evaluate(a.nrows); cols.evaluate(a.ncols);
	if (a.shape == Matrix_Shape::FULL) {
		double *vals = a.vals + rows.start * a.ncols + cols.start;
		cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans, rows.end - rows.start, cols.end - cols.start, alpha, vals, a.ncols, x.vals, 1, beta, y.vals, 1);
	} else {
		if (rows.start != cols.start || rows.end != cols.end) {
			eslog::error("cannot apply non-square slice to a packed matrix.\n");
		}
		double *vals = a.vals + rows.start * a.ncols + cols.start - (rows.start) * (rows.end - rows.start + 1) / 2;
		CBLAS_UPLO uplo = a.shape == Matrix_Shape::UPPER ? CBLAS_UPLO::CblasUpper : CBLAS_UPLO::CblasLower;
		cblas_dspmv(CBLAS_LAYOUT::CblasRowMajor, uplo, rows.end - rows.start, alpha, vals, x.vals, 1, beta, y.vals, 1);
	}
}

template <>
void applyT(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_Dense<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	if (a.submatrix[0].step != 1 || a.submatrix[1].step != 1) {
		eslog::error("slice is incompatible with apply.");
	}
	Slice rows = a.submatrix[0], cols = a.submatrix[1];
	rows.evaluate(a.nrows); cols.evaluate(a.ncols);
	if (a.shape == Matrix_Shape::FULL) {
		std::complex<double> *vals = a.vals + rows.start * a.ncols + cols.start;
		cblas_zgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans, a.nrows, a.ncols, &alpha, vals, a.ncols, x.vals, 1, &beta, y.vals, 1);
	} else {
		eslog::error("not implemented BLAS routine.\n");
	}
}

}
}

#endif
#endif
