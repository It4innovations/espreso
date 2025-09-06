
#include "esinfo/eslog.h"
#include "math/wrappers/math.blas.h"

#include <complex>

#ifdef HAVE_BLAS
#ifdef ESPRESO_USE_WRAPPER_DNBLAS_BLAS

#include "cblas.h"

#if defined(BLIS_INT_TYPE_SIZE) or defined(NVPL_BLAS_API)
using CBLAS_LAYOUT = CBLAS_ORDER;
#endif

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

template <typename T, typename I>
void apply(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x)
{
    if constexpr(std::is_same_v<T,double>) {
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
    else if constexpr(std::is_same_v<T,std::complex<double>>) {
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
    else {
        eslog::error("blas apply: unimplemented instantiation");
    }
}

template <typename T, typename I>
void applyT(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, const T &beta, const Vector_Dense<T, I> &x)
{
    if constexpr(std::is_same_v<T,double>) {
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
    else if constexpr(std::is_same_v<T,std::complex<double>>) {
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
    else {
        eslog::error("blas applyT: unimplemented instantiation");
    }
}

template <typename T, typename I>
void apply_hermitian(Vector_Dense<T, I> &y, const T &alpha, const Matrix_Dense<T, I> &a, char uplo, const T &beta, const Vector_Dense<T, I> &x)
{
    if(getSymmetry(a.type) != Matrix_Symmetry::HERMITIAN) eslog::error("blas apply_hermitian: matrix has to be hermitian\n");
    CBLAS_UPLO uplo_val;
    if(uplo == 'U') uplo_val = CblasUpper;
    else if(uplo == 'L') uplo_val = CblasLower;
    else eslog::error("blas apply_hermitian: invalid uplo\n");
    if(a.nrows != a.ncols) eslog::error("blas apply_hermitian: matrix is not square\n");
    if(x.size != a.ncols || y.size != a.nrows) eslog::error("blas apply_hermitian: incompatible vector size\n");
    if constexpr(std::is_same_v<T, float>)                cblas_ssymv(CblasRowMajor, uplo_val, a.nrows,  alpha, a.vals, a.get_ld(), x.vals, 1,  beta, y.vals, 1);
    if constexpr(std::is_same_v<T, double>)               cblas_dsymv(CblasRowMajor, uplo_val, a.nrows,  alpha, a.vals, a.get_ld(), x.vals, 1,  beta, y.vals, 1);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_chemv(CblasRowMajor, uplo_val, a.nrows, &alpha, a.vals, a.get_ld(), x.vals, 1, &beta, y.vals, 1);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zhemv(CblasRowMajor, uplo_val, a.nrows, &alpha, a.vals, a.get_ld(), x.vals, 1, &beta, y.vals, 1);
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

template<typename T>
void trsm(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & X, T alpha)
{
    if(A.nrows != A.ncols) eslog::error("system matrix has to be square\n");
    if(X.nrows != A.nrows) eslog::error("matrices are incompatible\n");

    auto layout = ((X.order == 'R') ? CblasRowMajor : CblasColMajor);
    auto trans_A = ((A.order == X.order) ? CblasNoTrans : CblasTrans);
    auto uplo_A = (((A.prop.uplo == 'U') == (A.order == X.order)) ? CblasUpper : CblasLower);
    auto diag_A = ((A.prop.diag == 'U') ? CblasUnit : CblasNonUnit);

    if constexpr(std::is_same_v<T, float>)                cblas_strsm(layout, CblasLeft, uplo_A, trans_A, diag_A, X.nrows, X.ncols,  alpha, A.vals, A.ld, X.vals, X.ld);
    if constexpr(std::is_same_v<T, double>)               cblas_dtrsm(layout, CblasLeft, uplo_A, trans_A, diag_A, X.nrows, X.ncols,  alpha, A.vals, A.ld, X.vals, X.ld);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_ctrsm(layout, CblasLeft, uplo_A, trans_A, diag_A, X.nrows, X.ncols, &alpha, A.vals, A.ld, X.vals, X.ld);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_ztrsm(layout, CblasLeft, uplo_A, trans_A, diag_A, X.nrows, X.ncols, &alpha, A.vals, A.ld, X.vals, X.ld);
}

template<typename T>
void gemm(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & B, MatrixDenseView_new<T> & C, T alpha, T beta)
{
    if(A.nrows != C.nrows || B.ncols != C.ncols || A.ncols != B.nrows) eslog::error("incompatible matrices");

    auto layout = ((C.order == 'R') ? CblasRowMajor : CblasColMajor);
    auto trans_A = ((A.order == C.order) ? CblasNoTrans : CblasTrans);
    auto trans_B = ((B.order == C.order) ? CblasNoTrans : CblasTrans);

    if constexpr(std::is_same_v<T, float>)                cblas_sgemm(layout, trans_A, trans_B, C.nrows, C.ncols, A.ncols,  alpha, A.vals, A.ld, B.vals, B.ld,  beta, C.vals, C.ld);
    if constexpr(std::is_same_v<T, double>)               cblas_dgemm(layout, trans_A, trans_B, C.nrows, C.ncols, A.ncols,  alpha, A.vals, A.ld, B.vals, B.ld,  beta, C.vals, C.ld);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_cgemm(layout, trans_A, trans_B, C.nrows, C.ncols, A.ncols, &alpha, A.vals, A.ld, B.vals, B.ld, &beta, C.vals, C.ld);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zgemm(layout, trans_A, trans_B, C.nrows, C.ncols, A.ncols, &alpha, A.vals, A.ld, B.vals, B.ld, &beta, C.vals, C.ld);
}

template<typename T>
void herk(MatrixDenseView_new<T> & A, MatrixDenseView_new<T> & C, herk_mode mode, utils::remove_complex_t<T> alpha, utils::remove_complex_t<T> beta)
{
    if(C.nrows != C.ncols) eslog::error("matrix C must be square\n");
    if(C.prop.uplo != 'U' && C.prop.uplo != 'L') eslog::error("C does not have set uplo\n");
    if(mode == herk_mode::AAh && A.nrows != C.nrows) eslog::error("incompatible matrix sizes\n");
    if(mode == herk_mode::AhA && A.ncols != C.ncols) eslog::error("incompatible matrix sizes\n");

    size_t n = C.nrows;
    size_t k = ((mode == herk_mode::AAh) ? A.ncols : A.nrows);

    auto layout = ((C.order == 'R') ? CblasRowMajor : CblasColMajor);
    auto uplo = ((C.prop.uplo == 'U') ? CblasUpper : CblasLower);
    auto trans = (((A.order == C.order) == (mode == herk_mode::AAh)) ? CblasNoTrans : CblasConjTrans);
    bool need_conj = (utils::is_complex<T>() && (A.order != C.order));

    if constexpr(utils::is_complex<T>()) if(need_conj) matrix_conj(C.vals, C.nrows, C.ncols, C.ld, C.order, C.prop.uplo);
    if constexpr(std::is_same_v<T, float>)                cblas_ssyrk(layout, uplo, trans, n, k, alpha, A.vals, A.ld, beta, C.vals, C.ld);
    if constexpr(std::is_same_v<T, double>)               cblas_dsyrk(layout, uplo, trans, n, k, alpha, A.vals, A.ld, beta, C.vals, C.ld);
    if constexpr(std::is_same_v<T, std::complex<float>>)  cblas_cherk(layout, uplo, trans, n, k, alpha, A.vals, A.ld, beta, C.vals, C.ld);
    if constexpr(std::is_same_v<T, std::complex<double>>) cblas_zherk(layout, uplo, trans, n, k, alpha, A.vals, A.ld, beta, C.vals, C.ld);
    if constexpr(utils::is_complex<T>()) if(need_conj) matrix_conj(C.vals, C.nrows, C.ncols, C.ld, C.order, C.prop.uplo);
}

template<typename T, bool conj>
static void transpose_internal(size_t src_nrows, size_t src_ncols, const T * src, size_t src_ld, T * dst, size_t dst_ld)
{
    // just a very basic tiled implementation. not tested.
    constexpr size_t tile_size = 64;
    size_t num_tiles_row = (src_nrows - 1) / tile_size + 1;
    size_t num_tiles_col = (src_ncols - 1) / tile_size + 1;
    for(size_t tile_row = 0; tile_row < num_tiles_row; tile_row++) {
        size_t r_start = tile_row * tile_size;
        size_t r_end = std::min(r_start + tile_size, src_nrows);
        size_t curr_tile_nrows = r_end - r_start;
        for(size_t tile_col = 0; tile_col < num_tiles_col; tile_col++) {
            size_t c_start = tile_col * tile_size;
            size_t c_end = std::min(c_start + tile_size, src_ncols);
            size_t curr_tile_ncols = c_end - c_start;
            const T * sub_src = src + r_start * src_ld + c_start;
            T * sub_dst = dst + c_start * dst_ld + r_start;
            for(size_t r = 0; r < curr_tile_nrows; r++) {
                for(size_t c = 0; c < curr_tile_ncols; c++) {
                    if constexpr(conj) {
                        sub_dst[c * dst_ld + r] = std::conj(sub_src[r * src_ld + c]);
                    }
                    else {
                        sub_dst[c * dst_ld + r] = sub_src[r * src_ld + c];
                    }
                }
            }
        }
    }
}

template<typename T>
void transpose(size_t src_nrows, size_t src_ncols, const T * src, size_t src_ld, T * dst, size_t dst_ld, char order, bool conj)
{
    if(order == 'C') {
        transpose<T>(src_ncols, src_nrows, src, src_ld, dst, dst_ld, 'R', conj);
        return;
    }

    if constexpr(utils::is_complex<T>()) {
        if(conj) {
            transpose_internal<T,true>(src_ncols, src_nrows, src, src_ld, dst, dst_ld);
        }
        else {
            transpose_internal<T,false>(src_ncols, src_nrows, src, src_ld, dst, dst_ld);
        }
    }
    else {
        transpose_internal<T,false>(src_ncols, src_nrows, src, src_ld, dst, dst_ld);
    }
}

template<typename T>
void transpose_inplace(size_t size, T * matrix, size_t ld, char order, bool conj)
{
    // just do it out-of-place and copy. todo better.

    constexpr size_t align = 64 / sizeof(T);
    size_t tmp_ld = ((size - 1) / align + 1) * align;

    T * tmp = new T[size * tmp_ld];

    transpose(size, size, matrix, ld, tmp, tmp_ld, order, conj);
    for(size_t i = 0; i < size; i++) {
        std::copy_n(tmp + i * tmp_ld, size, matrix + i * ld);
    }

    delete[] tmp;
}

template<typename T>
void hemm(size_t m, size_t n, T alpha, T * A, char order_A, size_t lda, char uplo_A, T * B, char order_B, size_t ldb, T beta, T * C, char order_C, size_t ldc)
{
    if(order_B != order_C) eslog::error("orders of B and C must match\n");

    if(order_A != order_C) {
        if constexpr(utils::is_complex<T>()) eslog::error("for complex matrices, all orders must match\n"); // todo do some copies and conjugations
        order_A = order_C;
        uplo_A = ((uplo_A == 'L') ? 'U' : 'L');
    }

    auto layout = ((order_C == 'R') ? CblasRowMajor : CblasColMajor);
    auto side = CblasLeft;
    auto uplo = ((uplo_A == 'U') ? CblasUpper : CblasLower);

    if constexpr(std::is_same_v<T,float>)                cblas_ssymm(layout, side, uplo, m, n,  alpha, A, lda, B, ldb,  beta, C, ldc);
    if constexpr(std::is_same_v<T,double>)               cblas_dsymm(layout, side, uplo, m, n,  alpha, A, lda, B, ldb,  beta, C, ldc);
    if constexpr(std::is_same_v<T,std::complex<float>>)  cblas_chemm(layout, side, uplo, m, n, &alpha, A, lda, B, ldb, &beta, C, ldc);
    if constexpr(std::is_same_v<T,std::complex<double>>) cblas_zhemm(layout, side, uplo, m, n, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

}
}
}

#include "math/wrappers/math.blas.inst.hpp"

#endif
#endif
