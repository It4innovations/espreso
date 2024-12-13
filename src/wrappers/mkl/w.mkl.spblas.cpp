
#include "math/math.h"
#include "math/primitives/matrix_csr.h"
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_MKL
#ifdef ESPRESO_USE_WRAPPER_SPBLAS_MKL

#include "mkl_spblas.h"

namespace espreso {

struct Matrix_SpBLAS_External_Representation {
    sparse_matrix_t inspector;
};

static void checkStatus(sparse_status_t status)
{
    switch (status) {
    case SPARSE_STATUS_SUCCESS: break;
    case SPARSE_STATUS_NOT_INITIALIZED: eslog::error("MKL SPARSE_STATUS_NOT_INITIALIZED\n"); break;
    case SPARSE_STATUS_ALLOC_FAILED: eslog::error("MKL SPARSE_STATUS_ALLOC_FAILED\n"); break;
    case SPARSE_STATUS_INVALID_VALUE: eslog::error("MKL SPARSE_STATUS_INVALID_VALUE\n"); break;
    case SPARSE_STATUS_EXECUTION_FAILED: eslog::error("MKL SPARSE_STATUS_EXECUTION_FAILED\n"); break;
    case SPARSE_STATUS_INTERNAL_ERROR: eslog::error("MKL SPARSE_STATUS_INTERNAL_ERROR\n"); break;
    case SPARSE_STATUS_NOT_SUPPORTED: eslog::error("MKL SPARSE_STATUS_NOT_SUPPORTED\n"); break;
}
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::~SpBLAS()
{
    if (_spblas) {
        checkStatus(mkl_sparse_destroy(_spblas->inspector));
    }
}

template <template <typename, typename, typename> class Matrix, typename T, typename I> void create(const typename SpBLAS<Matrix,T,I>::MatrixType *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed = false)
{
    eslog::error("SpBLAS wrapper is incompatible with T=%dB, I=%dB\n", sizeof(T), sizeof(I));
}

template <> void create<Matrix_CSR, float, int>(const Matrix_CSR<float, int> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed)
{
    if (matrix->nnz) {
        if (transposed) {
            switch (Indexing::CSR) {
            case 0: checkStatus(mkl_sparse_s_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            case 1: checkStatus(mkl_sparse_s_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            }
        } else {
            switch (Indexing::CSR) {
            case 0: checkStatus(mkl_sparse_s_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            case 1: checkStatus(mkl_sparse_s_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            }
        }
    }
}

template <> void create<Matrix_CSR, double, int>(const Matrix_CSR<double, int> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed)
{
    if (matrix->nnz) {
        if (transposed) {
            switch (Indexing::CSR) {
            case 0: checkStatus(mkl_sparse_d_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->ncols, matrix->nrows, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            case 1: checkStatus(mkl_sparse_d_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->ncols, matrix->nrows, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            }
        } else {
            switch (Indexing::CSR) {
            case 0: checkStatus(mkl_sparse_d_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            case 1: checkStatus(mkl_sparse_d_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            }
        }

    }
}

template <> void create<Matrix_CSR, std::complex<double>, int>(const Matrix_CSR<std::complex<double>, int> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed)
{
    if (matrix->nnz) {
        if (transposed) {
            switch (Indexing::CSR) {
            case 0: checkStatus(mkl_sparse_z_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            case 1: checkStatus(mkl_sparse_z_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            }
        } else {
            switch (Indexing::CSR) {
            case 0: checkStatus(mkl_sparse_z_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            case 1: checkStatus(mkl_sparse_z_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
            }
        }
    }
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS(MatrixType &a, bool trans)
: matrix{&a}, _spblas{}
{
    _spblas = new Matrix_SpBLAS_External_Representation();
    create<Matrix, T, I>(matrix, _spblas, trans);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(MatrixType &a, bool trans)
{
    matrix = &a;
    if (_spblas) {
        checkStatus(mkl_sparse_destroy(_spblas->inspector));
    }
    _spblas = new Matrix_SpBLAS_External_Representation();
    create<Matrix, T, I>(matrix, _spblas, trans);
}

static void setDescription(matrix_descr &descr, Matrix_Type type, Matrix_Shape shape)
{
    switch (type) {
    case Matrix_Type::UNSET_INVALID_NONE:
        eslog::error("Invalid/unset matrix type\n");
        break;
    case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
    case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:
    case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
    case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:
    case Matrix_Type::COMPLEX_SYMMETRIC:
        descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        break;
    case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:
    case Matrix_Type::REAL_NONSYMMETRIC:
    case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:
    case Matrix_Type::COMPLEX_NONSYMMETRIC:
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        break;
    }

    switch (shape) {
    case Matrix_Shape::FULL : descr.mode = SPARSE_FILL_MODE_FULL; break;
    case Matrix_Shape::UPPER: descr.mode = SPARSE_FILL_MODE_UPPER; break;
    case Matrix_Shape::LOWER: descr.mode = SPARSE_FILL_MODE_LOWER; break;
    }
}

template <>
void SpBLAS<Matrix_CSR, float, int>::apply(Vector_Dense<float> &y, const float &alpha, const float &beta, const Vector_Dense<float> &x)
{
    matrix_descr descr;
    setDescription(descr, matrix->type, matrix->shape);
    descr.diag = SPARSE_DIAG_NON_UNIT;
    checkStatus(mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{
    matrix_descr descr;
    setDescription(descr, matrix->type, matrix->shape);
    descr.diag = SPARSE_DIAG_NON_UNIT;
    checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void SpBLAS<Matrix_CSR, std::complex<double>, int>::apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
    matrix_descr descr;
    setDescription(descr, matrix->type, matrix->shape);
    descr.diag = SPARSE_DIAG_NON_UNIT;
    checkStatus(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Vector_Dense<T, I> &y, const T &alpha, const T &beta, const Vector_Dense<T, I> &x)
{
    eslog::error("SpBLAS wrapper is incompatible with T=%dB, I=%dB\n", sizeof(T), sizeof(I));
}

template <>
void SpBLAS<Matrix_CSR, float, int>::apply(Matrix_Dense<float> &y, const float &alpha, const float &beta, const Matrix_Dense<float> &x, bool trans)
{
    matrix_descr descr;
    setDescription(descr, matrix->type, matrix->shape);
    descr.diag = SPARSE_DIAG_NON_UNIT;
    checkStatus(mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void SpBLAS<Matrix_CSR, double, int>::apply(Matrix_Dense<double> &y, const double &alpha, const double &beta, const Matrix_Dense<double> &x, bool trans)
{
    matrix_descr descr;
    setDescription(descr, matrix->type, matrix->shape);
    descr.diag = SPARSE_DIAG_NON_UNIT;
    sparse_layout_t layout = trans ? SPARSE_LAYOUT_COLUMN_MAJOR : SPARSE_LAYOUT_ROW_MAJOR;
    checkStatus(mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, layout, x.vals, y.ncols, x.ncols, beta, y.vals, y.ncols));
}

template <>
void SpBLAS<Matrix_CSR, std::complex<double>, int>::apply(Matrix_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Matrix_Dense<std::complex<double> > &x, bool trans)
{
    matrix_descr descr;
    setDescription(descr, matrix->type, matrix->shape);
    descr.diag = SPARSE_DIAG_NON_UNIT;
    checkStatus(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Matrix_Dense<T, I> &y, const T &alpha, const T &beta, const Matrix_Dense<T, I> &x, bool trans)
{
    eslog::error("SpBLAS wrapper is incompatible with T=%dB, I=%dB\n", sizeof(T), sizeof(I));
}

}

#include "math/wrappers/math.spblas.inst.hpp"

#endif
#endif
