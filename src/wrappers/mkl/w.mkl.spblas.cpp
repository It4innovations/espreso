
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

namespace math {
namespace spblas {

class mkl_sparse_matrix_wrapper
{
public:
    sparse_matrix_t mat;
public:
    template<typename T, typename I>
    mkl_sparse_matrix_wrapper(MatrixCsxView_new<T,I> & A)
    {
        if(Indexing::CSR != 0) eslog::error("i dont support one-based csr indexing\n");
        if constexpr(std::is_same_v<T,float>)                checkStatus(mkl_sparse_s_create_csr(&mat, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
        if constexpr(std::is_same_v<T,double>)               checkStatus(mkl_sparse_d_create_csr(&mat, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
        if constexpr(std::is_same_v<T,std::complex<float>>)  checkStatus(mkl_sparse_c_create_csr(&mat, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
        if constexpr(std::is_same_v<T,std::complex<double>>) checkStatus(mkl_sparse_z_create_csr(&mat, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
    }
    ~mkl_sparse_matrix_wrapper()
    {
        checkStatus(mkl_sparse_destroy(mat));
    }
};

struct _handle_trsm
{
    std::unique_ptr<mkl_sparse_matrix_wrapper> A_mkl;
    sparse_operation_t op;
    matrix_descr A_descr;
    sparse_layout_t dense_layout;
    void * A_vals;
};

template<typename T, typename I>
void create_mkl_csr_matrix(sparse_matrix_t & A_mkl, MatrixCsxView_new<T,I> & A)
{
    if(Indexing::CSR != 0) eslog::error("i dont support one-based csr indexing\n");
    if constexpr(std::is_same_v<T,float>)                checkStatus(mkl_sparse_s_create_csr(&A_mkl, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
    if constexpr(std::is_same_v<T,double>)               checkStatus(mkl_sparse_d_create_csr(&A_mkl, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
    if constexpr(std::is_same_v<T,std::complex<float>>)  checkStatus(mkl_sparse_c_create_csr(&A_mkl, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
    if constexpr(std::is_same_v<T,std::complex<double>>) checkStatus(mkl_sparse_z_create_csr(&A_mkl, SPARSE_INDEX_BASE_ZERO, A.get_size_primary(), A.get_size_secdary(), A.ptrs, A.ptrs+1, A.idxs, A.vals));
}

template<typename T, typename I>
void trsm(MatrixCsxView_new<T,I> & A, MatrixDenseView_new<T> & X, MatrixDenseView_new<T> & Y, handle_trsm & handle, char stage)
{
    if(A.nrows != A.ncols) eslog::error("system matrix A has to be square\n");
    if(X.nrows != Y.nrows || X.ncols != Y.ncols) eslog::error("X and Y sizes dont match\n");
    if(X.order != Y.order) eslog::error("X and Y orders must match\n");
    if(A.nrows != X.nrows) eslog::error("incompatible matrix sizes\n");

    if(stage == 'A') { // All at once
        trsm<T,I>(A, X, Y, handle, 'p'); // no need to optimize for a single call
        trsm<T,I>(A, X, Y, handle, 'C');
    }
    if(stage == 'P' || stage == 'p') { // Preprocess
        if(handle.get() != nullptr) eslog::error("preprocess was already called\n");
        handle = std::make_shared<_handle_trsm>();
        handle->A_vals = A.vals;
        handle->A_mkl = std::make_unique<mkl_sparse_matrix_wrapper>(A);
        if(A.order == 'R') handle->op = SPARSE_OPERATION_NON_TRANSPOSE;
        if(A.order == 'C') handle->op = SPARSE_OPERATION_TRANSPOSE;
        handle->A_descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
        if((A.prop.uplo == 'U') == (A.order == 'R')) handle->A_descr.mode = SPARSE_FILL_MODE_UPPER;
        if((A.prop.uplo == 'L') == (A.order == 'R')) handle->A_descr.mode = SPARSE_FILL_MODE_LOWER;
        if(A.prop.diag == 'U') handle->A_descr.diag = SPARSE_DIAG_UNIT;
        if(A.prop.diag == 'N') handle->A_descr.diag = SPARSE_DIAG_NON_UNIT;
        if(X.order == 'R') handle->dense_layout = SPARSE_LAYOUT_ROW_MAJOR;
        if(X.order == 'C') handle->dense_layout = SPARSE_LAYOUT_COLUMN_MAJOR;
        if(stage == 'P') { // optimize only for capital P
            mkl_sparse_set_sm_hint(handle->A_mkl, handle->op, handle->A_descr, handle->dense_layout, X.ncols, 1000);
            checkStatus(mkl_sparse_optimize(handle->A_mkl));
        }
    }
    if(stage == 'C') { // Compute
        if(handle.get() == nullptr) eslog::error("preprocess has not been called\n");
        if(handle->A_vals != A.vals) eslog::error("matrix reallocation is not supported\n");
        if constexpr(std::is_same_v<T,float>)                checkStatus(mkl_sparse_s_trsm(handle->op, 1.0, handle->A_mkl, handle->A_descr, handle->dense_layout, X.vals, X.ncols, X.ld, Y.vals, Y.ld));
        if constexpr(std::is_same_v<T,double>)               checkStatus(mkl_sparse_d_trsm(handle->op, 1.0, handle->A_mkl, handle->A_descr, handle->dense_layout, X.vals, X.ncols, X.ld, Y.vals, Y.ld));
        if constexpr(std::is_same_v<T,std::complex<float>>)  checkStatus(mkl_sparse_c_trsm(handle->op, 1.0, handle->A_mkl, handle->A_descr, handle->dense_layout, X.vals, X.ncols, X.ld, Y.vals, Y.ld));
        if constexpr(std::is_same_v<T,std::complex<double>>) checkStatus(mkl_sparse_z_trsm(handle->op, 1.0, handle->A_mkl, handle->A_descr, handle->dense_layout, X.vals, X.ncols, X.ld, Y.vals, Y.ld));
    }
}

struct _handle_mm
{
    std::unique_ptr<mkl_sparse_matrix_wrapper> A_mkl;
    sparse_operation_t op;
    matrix_descr A_descr;
    sparse_layout_t dense_layout;
    void * A_vals;
};

template<typename T, typename I>
void mm(MatrixCsxView_new<T,I> & A, MatrixDenseView_new<T> & B, MatrixDenseView_new<T> & C, T alpha, T beta, handle_mm & handle, char stage)
{
    if(A.nrows != C.nrows || B.ncols != C.ncols || A.ncols != B.nrows) eslog::error("incompatible matrices\n");
    if(B.order != C.order) eslog::error("B and C order must match\n");

    if(stage == 'A') { // All at once
        mm<T,I>(A, B, C, alpha, beta, handle, 'p'); // no need to optimize for a single call
        mm<T,I>(A, B, C, alpha, beta, handle, 'C');
    }
    if(stage == 'P' || stage == 'p') { // Preprocess
        if(handle.get() != nullptr) eslog::error("preprocess was already called\n");
        handle = std::make_shared<_handle_mm>();
        handle->A_vals = A.vals;
        handle->A_mkl = std::make_unique<mkl_sparse_matrix_wrapper>(A);
        handle->op = ((A.order == 'R') ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE);
        handle->A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        handle->dense_layout = ((B.order == 'R') ? SPARSE_LAYOUT_ROW_MAJOR : SPARSE_LAYOUT_COLUMN_MAJOR);
        if(stage == 'P') { // optimize only for capical P
            mkl_sparse_set_mm_hint(handle->A_mkl, handle->op, handle->A_descr, handle->dense_layout, B.ncols, 1000);
            checkStatus(mkl_sparse_optimize(handle->A_mkl));
        }
    }
    if(stage == 'C') { // Compute
        if(handle.get() == nullptr) eslog::error("preprocess has not been called\n");
        if(handle->A_vals != A.vals) eslog::error("matrix reallocation is not supported\n");
        if constexpr(std::is_same_v<T,float>)                checkStatus(mkl_sparse_s_mm(handle->op, alpha, handle->A_mkl, handle->A_descr, handle->dense_layout, B.vals, B.ncols, B.ld, beta, C.vals, C.ld));
        if constexpr(std::is_same_v<T,double>)               checkStatus(mkl_sparse_d_mm(handle->op, alpha, handle->A_mkl, handle->A_descr, handle->dense_layout, B.vals, B.ncols, B.ld, beta, C.vals, C.ld));
        if constexpr(std::is_same_v<T,std::complex<float>>)  checkStatus(mkl_sparse_c_mm(handle->op, alpha, handle->A_mkl, handle->A_descr, handle->dense_layout, B.vals, B.ncols, B.ld, beta, C.vals, C.ld));
        if constexpr(std::is_same_v<T,std::complex<double>>) checkStatus(mkl_sparse_z_mm(handle->op, alpha, handle->A_mkl, handle->A_descr, handle->dense_layout, B.vals, B.ncols, B.ld, beta, C.vals, C.ld));
    }
}

}
}

}

#include "math/wrappers/math.spblas.inst.hpp"

#endif
#endif
