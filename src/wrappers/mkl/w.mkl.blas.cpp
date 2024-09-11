
#include "esinfo/eslog.h"
#include "math/wrappers/math.blas.h"

#include <complex>

#ifdef HAVE_MKL
#include "mkl_blas.h"
#include "mkl_cblas.h"

namespace espreso {
namespace math {
namespace blas {

template <typename T>
void copy(const int size, T *x, const int incX, const T *y, const int incY)
{
    if constexpr(std::is_same_v<T, float>)                cblas_scopy(size, y, incY, x, incX);
    if constexpr(std::is_same_v<T, double>)               cblas_dcopy(size, y, incY, x, incX);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_ccopy(size, y, incY, x, incX);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zcopy(size, y, incY, x, incX);
}

template <typename T>
void scale(const int size, const T &alpha, T *x, const int incX)
{
    if constexpr(std::is_same_v<T, float>)                cblas_sscal(size,  alpha, x, incX);
    if constexpr(std::is_same_v<T, double>)               cblas_dscal(size,  alpha, x, incX);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_cscal(size, &alpha, x, incX);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zscal(size, &alpha, x, incX);
}

template <typename T>
void add(const int size, T *x, const int incX, const T &alpha, const T *y, const int incY)
{
    if constexpr(std::is_same_v<T, float>)                cblas_saxpy(size,  alpha, y, incY, x, incX);
    if constexpr(std::is_same_v<T, double>)               cblas_daxpy(size,  alpha, y, incY, x, incX);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_caxpy(size, &alpha, y, incY, x, incX);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zaxpy(size, &alpha, y, incY, x, incX);
}

template <typename T>
T dot(const int size, const T *x, const int incX, const T *y, const int incY)
{
    if constexpr(std::is_same_v<T, float>)                return cblas_sdot(size, x, incX, y, incY);
    if constexpr(std::is_same_v<T, double>)               return cblas_ddot(size, x, incX, y, incY);
    if constexpr(std::is_same_v<T, std::complex<float>>)  eslog::error("not implemented BLAS routine.\n");
    if constexpr(std::is_same_v<T, std::complex<double>>) eslog::error("not implemented BLAS routine.\n");
//    return cblas_cdotu_sub(size, x, incX, y, incY, &dot);
}

template <typename T>
utils::remove_complex_t<T> norm(const int size, const T *x, const int incX)
{
    if constexpr(std::is_same_v<T, float>)                return cblas_snrm2(size, x, incX);
    if constexpr(std::is_same_v<T, double>)               return cblas_dnrm2(size, x, incX);
    if constexpr(std::is_same_v<T, std::complex<float>>)  return cblas_scnrm2(size, x, incX);
    if constexpr(std::is_same_v<T, std::complex<double>>) return cblas_dznrm2(size, x, incX);
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

template <typename T, typename I>
void apply(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("blas apply: unimplemented instantiation");
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

template <typename T, typename I>
void applyT(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("blas applyT: unimplemented instantiation");
}

template <typename T, typename I>
void apply_hermitian(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x)
{
    if(getSymmetry(a.type) != Matrix_Symmetry::HERMITIAN) eslog::error("blas apply_hermitian: matrix has to be hermitian\n");
    CBLAS_UPLO uplo;
    if(a.shape == Matrix_Shape::UPPER) uplo = CblasUpper;
    else if(a.shape == Matrix_Shape::LOWER) uplo = CblasLower;
    else eslog::error("blas apply_hermitian: invalid matrix shape\n");
    if(a.nrows != a.ncols) eslog::error("blas apply_hermitian: matrix is not square\n");
    if(x.size != a.ncols || y.size != a.nrows) eslog::error("blas apply_hermitian: incompatible vector size\n");
    if constexpr(std::is_same_v<T, float>)                cblas_ssymv(CblasRowMajor, uplo, a.nrows,  alpha, a.vals, a.get_ld(), x.vals, 1,  beta, y.vals, 1);
    if constexpr(std::is_same_v<T, double>)               cblas_dsymv(CblasRowMajor, uplo, a.nrows,  alpha, a.vals, a.get_ld(), x.vals, 1,  beta, y.vals, 1);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_chemv(CblasRowMajor, uplo, a.nrows, &alpha, a.vals, a.get_ld(), x.vals, 1, &beta, y.vals, 1);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zhemv(CblasRowMajor, uplo, a.nrows, &alpha, a.vals, a.get_ld(), x.vals, 1, &beta, y.vals, 1);
}

template <typename T, typename I>
void AAt(const Matrix_Dense<T, I> &A, Matrix_Dense<T, I> &AAt, bool trans)
{
    I size = trans ? A.ncols : A.nrows;
    if (AAt.nrows != AAt.ncols) eslog::error("invalid AAt dimension.\n");
    if (size != AAt.ncols) eslog::error("invalid AAt dimension.\n");
    T one = 1;
    T zero = 0;
    if constexpr(std::is_same_v<T, float>)                cblas_ssyrk(CblasRowMajor, CblasUpper, trans ? CblasTrans : CblasNoTrans, AAt.nrows, trans ? A.nrows : A.ncols,  one, A.vals, A.ncols,  zero, AAt.vals, AAt.ncols);
    if constexpr(std::is_same_v<T, double>)               cblas_dsyrk(CblasRowMajor, CblasUpper, trans ? CblasTrans : CblasNoTrans, AAt.nrows, trans ? A.nrows : A.ncols,  one, A.vals, A.ncols,  zero, AAt.vals, AAt.ncols);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_csyrk(CblasRowMajor, CblasUpper, trans ? CblasTrans : CblasNoTrans, AAt.nrows, trans ? A.nrows : A.ncols, &one, A.vals, A.ncols, &zero, AAt.vals, AAt.ncols);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zsyrk(CblasRowMajor, CblasUpper, trans ? CblasTrans : CblasNoTrans, AAt.nrows, trans ? A.nrows : A.ncols, &one, A.vals, A.ncols, &zero, AAt.vals, AAt.ncols);
}

template <typename T, typename I>
void multiply(T alpha, const Matrix_Dense<T, I> &A, const Matrix_Dense<T, I> &B, T beta, Matrix_Dense<T, I> &C, bool transA, bool transB)
{
    I rows = transA ? A.ncols : A.nrows, cols = transB ? B.nrows : B.ncols;
    if (C.nrows != rows || C.ncols != cols) eslog::error("invalid dimension.\n");
    if (A.shape != Matrix_Shape::FULL) eslog::error("invalid shape.\n");
    if (B.shape != Matrix_Shape::FULL) eslog::error("invalid shape.\n");
    if constexpr(std::is_same_v<T, float>)                cblas_sgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, transB ? CblasTrans : CblasNoTrans, C.nrows, C.ncols, transA ? A.nrows : A.ncols,  alpha, A.vals, A.ncols, B.vals, B.ncols,  beta, C.vals, C.ncols);
    if constexpr(std::is_same_v<T, double>)               cblas_dgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, transB ? CblasTrans : CblasNoTrans, C.nrows, C.ncols, transA ? A.nrows : A.ncols,  alpha, A.vals, A.ncols, B.vals, B.ncols,  beta, C.vals, C.ncols);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_cgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, transB ? CblasTrans : CblasNoTrans, C.nrows, C.ncols, transA ? A.nrows : A.ncols, &alpha, A.vals, A.ncols, B.vals, B.ncols, &beta, C.vals, C.ncols);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, transB ? CblasTrans : CblasNoTrans, C.nrows, C.ncols, transA ? A.nrows : A.ncols, &alpha, A.vals, A.ncols, B.vals, B.ncols, &beta, C.vals, C.ncols);
}

template <typename T, typename I>
void multiply(T alpha, const Matrix_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, Vector_Dense<T, I> &C, bool transA)
{
    I rows = transA ? A.ncols : A.nrows;
    if (C.size != rows) eslog::error("invalid dimension.\n");
    if (A.shape != Matrix_Shape::FULL) eslog::error("invalid shape.\n");
    if constexpr(std::is_same_v<T, float>)                cblas_sgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, CblasNoTrans, C.size, 1, transA ? A.nrows : A.ncols,  alpha, A.vals, A.ncols, B.vals, 1,  beta, C.vals, 1);
    if constexpr(std::is_same_v<T, double>)               cblas_dgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, CblasNoTrans, C.size, 1, transA ? A.nrows : A.ncols,  alpha, A.vals, A.ncols, B.vals, 1,  beta, C.vals, 1);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_cgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, CblasNoTrans, C.size, 1, transA ? A.nrows : A.ncols, &alpha, A.vals, A.ncols, B.vals, 1, &beta, C.vals, 1);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, CblasNoTrans, C.size, 1, transA ? A.nrows : A.ncols, &alpha, A.vals, A.ncols, B.vals, 1, &beta, C.vals, 1);
}

template <typename T, typename I>
void multiply(T alpha, const Vector_Dense<T, I> &A, const Vector_Dense<T, I> &B, T beta, T &out)
{
    if (A.size != B.size) eslog::error("invalid dimension.\n");
    if constexpr(std::is_same_v<T, float>)                cblas_sgemv(CblasRowMajor, CblasNoTrans, 1, A.size,  alpha, A.vals, A.size, B.vals, 1,  beta, &out, 1);
    if constexpr(std::is_same_v<T, double>)               cblas_dgemv(CblasRowMajor, CblasNoTrans, 1, A.size,  alpha, A.vals, A.size, B.vals, 1,  beta, &out, 1);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_cgemv(CblasRowMajor, CblasNoTrans, 1, A.size, &alpha, A.vals, A.size, B.vals, 1, &beta, &out, 1);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zgemv(CblasRowMajor, CblasNoTrans, 1, A.size, &alpha, A.vals, A.size, B.vals, 1, &beta, &out, 1);
}

}
}
}

#include "math/wrappers/math.blas.inst.hpp"

#endif



