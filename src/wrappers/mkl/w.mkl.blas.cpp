
#include "esinfo/eslog.h"
#include "math/wrappers/math.blas.h"

#include <complex>

#ifdef HAVE_MKL
#include "mkl_blas.h"
#include "mkl_cblas.h"

namespace espreso {
namespace math {
namespace blas {

template <>
void copy(const int size, double *x, const int incX, const double *y, const int incY)
{
    cblas_dcopy(size, y, incY, x, incX);
}

template <>
void copy(const int size, std::complex<double> *x, const int incX, const std::complex<double> *y, const int incY)
{
    cblas_zcopy(size, y, incY, x, incX);
}

template <>
void scale(const int size, const float &alpha, float *x, const int incX)
{
    cblas_sscal(size, alpha, x, incX);
}

template <>
void scale(const int size, const double &alpha, double *x, const int incX)
{
    cblas_dscal(size, alpha, x, incX);
}

template <>
void scale(const int size, const std::complex<double> &alpha, std::complex<double> *x, const int incX)
{
    cblas_zscal(size, &alpha, x, incX);
}

template <>
void add(const int size, float *x, const int incX, const float &alpha, const float *y, const int incY)
{
    cblas_saxpy(size, alpha, y, incY, x, incX);
}

template <>
void add(const int size, double *x, const int incX, const double &alpha, const double *y, const int incY)
{
    cblas_daxpy(size, alpha, y, incY, x, incX);
}

template <>
void add(const int size, std::complex<float> *x, const int incX, const std::complex<float> &alpha, const std::complex<float> *y, const int incY)
{
    cblas_caxpy(size, &alpha, y, incY, x, incX);
}

template <>
void add(const int size, std::complex<double> *x, const int incX, const std::complex<double> &alpha, const std::complex<double> *y, const int incY)
{
    cblas_zaxpy(size, &alpha, y, incY, x, incX);
}

template <>
double dot(const int size, const double *x, const int incX, const double *y, const int incY)
{
    return cblas_ddot(size, x, incX, y, incY);
}

template <>
std::complex<double> dot(const int size, const std::complex<double> *x, const int incX, const std::complex<double> *y, const int incY)
{
    std::complex<double> dot = 0;
//    return cblas_cdotu_sub(size, x, incX, y, incY, &dot);
    return dot;
}

template <>
float norm(const int size, const float *x, const int incX)
{
    return cblas_snrm2(size, x, incX);
}

template <>
double norm(const int size, const double *x, const int incX)
{
    return cblas_dnrm2(size, x, incX);
}

template <>
float norm(const int size, const std::complex<float> *x, const int incX)
{
    return cblas_scnrm2(size, x, incX);
}

template <>
double norm(const int size, const std::complex<double> *x, const int incX)
{
    return cblas_dznrm2(size, x, incX);
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

template <>
void AAt(const Matrix_Dense<double> &A, Matrix_Dense<double> &AAt, bool trans)
{
    int size = trans ? A.ncols : A.nrows;
    if (AAt.nrows != AAt.ncols) eslog::error("invalid AAt dimension.\n");
    if (size != AAt.ncols) eslog::error("invalid AAt dimension.\n");
    cblas_dsyrk(CblasRowMajor, CblasUpper, trans ? CblasTrans : CblasNoTrans, AAt.nrows, trans ? A.nrows : A.ncols, 1, A.vals, A.ncols, 0, AAt.vals, AAt.ncols);
}

template <>
void multiply(double alpha, const Matrix_Dense<double> &A, const Matrix_Dense<double> &B, double beta, Matrix_Dense<double> &C, bool transA, bool transB)
{
    int rows = transA ? A.ncols : A.nrows, cols = transB ? B.nrows : B.ncols;
    if (C.nrows != rows || C.ncols != cols) eslog::error("invalid dimension.\n");
    cblas_dgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, transB ? CblasTrans : CblasNoTrans, C.nrows, C.ncols, transA ? A.nrows : A.ncols, alpha, A.vals, A.ncols, B.vals, B.ncols, beta, C.vals, C.ncols);
}

template <>
void multiply(double alpha, const Matrix_Dense<double> &A, const Vector_Dense<double> &B, double beta, Vector_Dense<double> &C, bool transA)
{
    int rows = transA ? A.ncols : A.nrows;
    if (C.size != rows) eslog::error("invalid dimension.\n");
    cblas_dgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, CblasNoTrans, C.size, 1, transA ? A.nrows : A.ncols, alpha, A.vals, A.ncols, B.vals, 1, beta, C.vals, 1);
}

template <>
void multiply(double alpha, const Vector_Dense<double> &A, const Vector_Dense<double> &B, double beta, double &out)
{
    if (A.size != B.size) eslog::error("invalid dimension.\n");
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 1, A.size, alpha, A.vals, A.size, B.vals, 1, beta, &out, 1);
}

template <>
void multiply(std::complex<double> alpha, const Vector_Dense<std::complex<double> > &A, const Vector_Dense<std::complex<double> > &B, std::complex<double> beta, std::complex<double> &out)
{
    if (A.size != B.size) eslog::error("invalid dimension.\n");
    cblas_zgemv(CblasRowMajor, CblasNoTrans, 1, A.size, &alpha, A.vals, A.size, B.vals, 1, &beta, &out, 1);
}

}
}
}

#endif



