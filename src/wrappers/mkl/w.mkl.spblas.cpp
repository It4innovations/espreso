
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_MKL
#include "mkl_spblas.h"

namespace espreso {

template struct SpBLAS<float, Matrix_CSR>;
template struct SpBLAS<double, Matrix_CSR>;
template struct SpBLAS<std::complex<double>, Matrix_CSR>;

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

template <typename T, template <typename> class Matrix>
SpBLAS<T, Matrix>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <typename T, template <typename> class Matrix>
SpBLAS<T, Matrix>::~SpBLAS()
{
	if (_spblas) {
		checkStatus(mkl_sparse_destroy(_spblas->inspector));
	}
}

template <typename T, template <typename> class Matrix> void create(const Matrix<T> *matrix, Matrix_SpBLAS_External_Representation *_spblas);

template <> void create<float, Matrix_CSR>(const Matrix_CSR<float> *matrix, Matrix_SpBLAS_External_Representation *_spblas)
{
	if (matrix->nnz) {
		switch (Indexing::CSR) {
		case 0: checkStatus(mkl_sparse_s_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
		case 1: checkStatus(mkl_sparse_s_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
		}
	}
}

template <> void create<double, Matrix_CSR>(const Matrix_CSR<double> *matrix, Matrix_SpBLAS_External_Representation *_spblas)
{
	if (matrix->nnz) {
		switch (Indexing::CSR) {
		case 0: checkStatus(mkl_sparse_d_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
		case 1: checkStatus(mkl_sparse_d_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
		}
	}
}

template <> void create<std::complex<double>, Matrix_CSR>(const Matrix_CSR<std::complex<double> > *matrix, Matrix_SpBLAS_External_Representation *_spblas)
{
	if (matrix->nnz) {
		switch (Indexing::CSR) {
		case 0: checkStatus(mkl_sparse_z_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
		case 1: checkStatus(mkl_sparse_z_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
		}
	}
}

template <typename T, template <typename> class Matrix>
SpBLAS<T, Matrix>::SpBLAS(const Matrix<T> &a)
: matrix{&a}, _spblas{}
{
	_spblas = new Matrix_SpBLAS_External_Representation();
	create(matrix, _spblas);
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::insert(const Matrix<T> &a)
{
	matrix = &a;
	if (_spblas) {
		checkStatus(mkl_sparse_destroy(_spblas->inspector));
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	create(matrix, _spblas);
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::insertTransposed(const Matrix<T> &a)
{
	eslog::error("MKL SpBLAS: call empty function.\n");
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::extractUpper(Matrix<T> &a)
{
	eslog::error("MKL SpBLAS: call empty function.\n");
}


static void setDescription(matrix_descr &descr, Matrix_Type type, Matrix_Shape shape)
{
	switch (type) {
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
void SpBLAS<float, Matrix_CSR>::apply(Vector_Dense<float> &y, const float &alpha, const float &beta, const Vector_Dense<float> &x)
{
	matrix_descr descr;
	setDescription(descr, matrix->type, matrix->shape);
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void SpBLAS<double, Matrix_CSR>::apply(Vector_Dense<double> &y, const double &alpha, const double &beta, const Vector_Dense<double> &x)
{
	matrix_descr descr;
	setDescription(descr, matrix->type, matrix->shape);
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void SpBLAS<std::complex<double>, Matrix_CSR>::apply(Vector_Dense<std::complex<double>> &y, const std::complex<double> &alpha, const std::complex<double> &beta, const Vector_Dense<std::complex<double>> &x)
{
	matrix_descr descr;
		setDescription(descr, matrix->type, matrix->shape);
		descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, _spblas->inspector, descr, x.vals, beta, y.vals));
}

template <typename T, template <typename> class Matrix>
void SpBLAS<T, Matrix>::transposeTo(SpBLAS<T, Matrix> &A)
{
	eslog::error("MKL SpBLAS: call empty function.\n");
}

template <>
void SpBLAS<float, Matrix_CSR>::multiply(SpBLAS<float, Matrix_CSR> &A, SpBLAS<float, Matrix_CSR> &B)
{
	eslog::error("MKL SpBLAS: call empty function.\n");
}

template <>
void SpBLAS<double, Matrix_CSR>::multiply(SpBLAS<double, Matrix_CSR> &A, SpBLAS<double, Matrix_CSR> &B)
{
	eslog::error("MKL SpBLAS: call empty function.\n");
}

}

#endif
