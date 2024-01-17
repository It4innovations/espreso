
#include "math/math.h"
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_MKL
#include "mkl_spblas.h"

namespace espreso {

template struct SpBLAS<Matrix_CSR, float, int>;
template struct SpBLAS<Matrix_CSR, double, int>;
template struct SpBLAS<Matrix_CSR, std::complex<double>, int>;

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

template <template <typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <template <typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::~SpBLAS()
{
	if (_spblas) {
		checkStatus(mkl_sparse_destroy(_spblas->inspector));
	}
}

template <template <typename, typename> class Matrix, typename T, typename I> void create(const Matrix<T, I> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed = false);
template <template <typename, typename> class Matrix, typename T, typename I> void extract(Matrix<T, I> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed = false);

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

template <> void extract<Matrix_CSR, float, int>(Matrix_CSR<float, int> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed)
{
//	if (matrix->nnz) {
//		if (transposed) {
//			switch (Indexing::CSR) {
//			case 0: checkStatus(mkl_sparse_s_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			case 1: checkStatus(mkl_sparse_s_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			}
//		} else {
//			switch (Indexing::CSR) {
//			case 0: checkStatus(mkl_sparse_s_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			case 1: checkStatus(mkl_sparse_s_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			}
//		}
//	}
}

template <> void extract<Matrix_CSR, double, int>(Matrix_CSR<double, int> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed)
{
	sparse_index_base_t indexing;
	esint *rowend;
	if (transposed) {
		checkStatus(mkl_sparse_d_export_csc(_spblas->inspector, &indexing, &matrix->nrows, &matrix->ncols, &matrix->rows, &rowend, &matrix->cols, &matrix->vals));
	} else {

		checkStatus(mkl_sparse_d_export_csr(_spblas->inspector, &indexing, &matrix->nrows, &matrix->ncols, &matrix->rows, &rowend, &matrix->cols, &matrix->vals));
		// TODO: check indexing
	}
	matrix->nnz = matrix->rows[matrix->nrows - 1 + matrix->rows[0]] - matrix->rows[0];
}

template <> void extract<Matrix_CSR, std::complex<double>, int>(Matrix_CSR<std::complex<double>, int> *matrix, Matrix_SpBLAS_External_Representation *_spblas, bool transposed)
{
//	if (matrix->nnz) {
//		if (transposed) {
//			switch (Indexing::CSR) {
//			case 0: checkStatus(mkl_sparse_z_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			case 1: checkStatus(mkl_sparse_z_create_csr(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			}
//		} else {
//			switch (Indexing::CSR) {
//			case 0: checkStatus(mkl_sparse_z_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ZERO, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			case 1: checkStatus(mkl_sparse_z_create_csc(&_spblas->inspector, SPARSE_INDEX_BASE_ONE, matrix->nrows, matrix->ncols, matrix->rows, matrix->rows + 1, matrix->cols, matrix->vals)); break;
//			}
//		}
//	}
}

template <template <typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS(Matrix<T, I> &a)
: matrix{&a}, _spblas{}
{
	_spblas = new Matrix_SpBLAS_External_Representation();
	create(matrix, _spblas, false);
}

template <template <typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(Matrix<T, I> &a)
{
	matrix = &a;
	if (_spblas) {
		checkStatus(mkl_sparse_destroy(_spblas->inspector));
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	create(matrix, _spblas, false);
}

template <>
void SpBLAS<Matrix_CSR, double, int>::insert(Matrix_Dense<double, int> &a, double threshold)
{
	matrix = new Matrix_CSR<double>();

	esint nnz = 0;
	for (esint r = 0; r < a.nrows; r++) {
		for (esint c = 0; c < a.ncols; c++) {
			if (a.vals[r * a.ncols + c] < -threshold || threshold < a.vals[r * a.ncols + c]) {
				++nnz;
			}
		}
	}
	matrix->resize(a.nrows, a.ncols, nnz);

	matrix->rows[0] = Indexing::CSR;
	for (esint r = 0, i = 0; r < a.nrows; r++) {
		for (esint c = 0; c < a.ncols; c++) {
			if (a.vals[r * a.ncols + c] < -threshold || threshold < a.vals[r * a.ncols + c]) {
				matrix->cols[i] = c + Indexing::CSR;
				matrix->vals[i] = a.vals[r * a.ncols + c];
				++i;
			}
		}
		matrix->rows[r + 1] = i + Indexing::CSR;
	}
	insert(*matrix);
}

template <template <typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insertTransposed(Matrix<T, I> &a)
{
	matrix = &a;
	if (_spblas) {
		checkStatus(mkl_sparse_destroy(_spblas->inspector));
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	create(matrix, _spblas, true);
	matrix = nullptr; // since it is not transposed
}

template <template <typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::extractUpper(Matrix<T, I> &a)
{
//	extract(&a, _spblas, false);
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

template <template <typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::transposeTo(SpBLAS<Matrix, T, I> &A)
{
	eslog::error("MKL SpBLAS: call empty function.\n");
}

template <template <typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::convertTo(Matrix_Dense<T, I> &out)
{
	// is there a better functions?
	out.resize(matrix->nrows, matrix->ncols);
	math::set(out, T{0.});

	for (esint r = 0; r < matrix->nrows; r++) {
		for (esint c = matrix->rows[r]; c < matrix->rows[r + 1]; ++c) {
			out.vals[r * matrix->ncols + matrix->cols[c - Indexing::CSR] - Indexing::CSR] = matrix->vals[c - Indexing::CSR];
		}
	}
}

template <>
void SpBLAS<Matrix_CSR, double, int>::multiply(SpBLAS<Matrix_CSR, double, int> &A, SpBLAS<Matrix_CSR, double, int> &B)
{
	if (_spblas) {
		checkStatus(mkl_sparse_destroy(_spblas->inspector));
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	checkStatus(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A._spblas->inspector, B._spblas->inspector, &_spblas->inspector));
	checkStatus(mkl_sparse_order(_spblas->inspector));
	matrix = new Matrix_CSR<double>();
	extract(matrix, _spblas, false);
}

template <>
void SpBLAS<Matrix_CSR, double, int>::multiply(SpBLAS<Matrix_CSR, double, int> &A, Matrix_Dense<double> &B)
{

}

template <>
void SpBLAS<Matrix_CSR, double, int>::AAt(SpBLAS<Matrix_CSR, double, int> &A)
{
	if (_spblas) {
		checkStatus(mkl_sparse_destroy(_spblas->inspector));
	}
	_spblas = new Matrix_SpBLAS_External_Representation();
	checkStatus(mkl_sparse_syrk(SPARSE_OPERATION_NON_TRANSPOSE, A._spblas->inspector, &_spblas->inspector));
	matrix = new Matrix_CSR<double>();
	extract(matrix, _spblas, false);
	matrix->shape = Matrix_Shape::UPPER;
	matrix->type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
}

template <>
void SpBLAS<Matrix_CSR, double, int>::solveRowMayor(Matrix_Dense<double> &rhs, Matrix_Dense<double> &solution)
{
	solution.resize(rhs.nrows, rhs.ncols);
	matrix_descr descr;
	setDescription(descr, matrix->type, matrix->shape);
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_trsm(SPARSE_OPERATION_NON_TRANSPOSE, 1., _spblas->inspector, descr, SPARSE_LAYOUT_ROW_MAJOR, rhs.vals, rhs.ncols, rhs.nrows, solution.vals, solution.ncols));
}

template <>
void SpBLAS<Matrix_CSR, double, int>::solveColMayor(Matrix_Dense<double> &rhs, Matrix_Dense<double> &solution)
{
	solution.resize(rhs.nrows, rhs.ncols);
	matrix_descr descr;
	setDescription(descr, matrix->type, matrix->shape);
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_trsm(SPARSE_OPERATION_NON_TRANSPOSE, 1., _spblas->inspector, descr, SPARSE_LAYOUT_COLUMN_MAJOR, rhs.vals, rhs.nrows, rhs.ncols, solution.vals, solution.nrows));
}

}

#endif
