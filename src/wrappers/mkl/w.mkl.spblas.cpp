
#include "math/wrappers/math.spblas.h"
#include "esinfo/eslog.h"

#include <complex>

#ifdef HAVE_MKL
#ifdef USE_SPBLAS_MKL
#include "mkl_spblas.h"

namespace espreso {

struct Matrix_CSR_External_Representation {
	sparse_matrix_t inspector;
};

struct Matrix_IJV_External_Representation {
	sparse_matrix_t inspector;
};

namespace math {

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

template <>
void commit(Matrix_Dense<double> &x)
{

}

template <>
void commit(Matrix_Dense<std::complex<double> > &x)
{

}

template <>
void commit(Matrix_CSR<double> &x)
{
	x._spblas = new Matrix_CSR_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_d_create_csr(&x._spblas->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.rows, x.rows + 1, x.cols, x.vals));
	}
}

template <>
void commit(Matrix_CSR<std::complex<double> > &x)
{
	x._spblas = new Matrix_CSR_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_z_create_csr(&x._spblas->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.rows, x.rows + 1, x.cols, x.vals));
	}
}

template <>
void commit(Matrix_IJV<double> &x)
{
	x._spblas = new Matrix_IJV_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_d_create_coo(&x._spblas->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.nnz, x.rows, x.cols, x.vals));
	}
}

template <>
void commit(Matrix_IJV<std::complex<double> > &x)
{
	x._spblas = new Matrix_IJV_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_z_create_coo(&x._spblas->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.nnz, x.rows, x.cols, x.vals));
	}
}

template <>
void free(Matrix_Dense<double> &x)
{

}

template <>
void free(Matrix_Dense<std::complex<double> > &x)
{

}

template <>
void free(Matrix_CSR<double> &x)
{
	if (x._spblas) {
		checkStatus(mkl_sparse_destroy(x._spblas->inspector));
		x._spblas = nullptr;
	}
}

template <>
void free(Matrix_CSR<std::complex<double> > &x)
{
	if (x._spblas) {
		checkStatus(mkl_sparse_destroy(x._spblas->inspector));
		x._spblas = nullptr;
	}
}

template <>
void free(Matrix_IJV<double> &x)
{
	if (x._spblas) {
		checkStatus(mkl_sparse_destroy(x._spblas->inspector));
		x._spblas = nullptr;
	}
}

template <>
void free(Matrix_IJV<std::complex<double> > &x)
{
	if (x._spblas) {
		checkStatus(mkl_sparse_destroy(x._spblas->inspector));
		x._spblas = nullptr;
	}
}


template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_CSR<double> &a, const double &beta, const Vector_Dense<double> &x)
{
	matrix_descr descr;
	switch (a.type) {
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

	switch (a.shape) {
	case Matrix_Shape::FULL : descr.mode = SPARSE_FILL_MODE_FULL; break;
	case Matrix_Shape::UPPER: descr.mode = SPARSE_FILL_MODE_UPPER; break;
	case Matrix_Shape::LOWER: descr.mode = SPARSE_FILL_MODE_LOWER; break;
	}
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_CSR<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	matrix_descr descr;
	switch (a.type) {
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

	switch (a.shape) {
	case Matrix_Shape::FULL : descr.mode = SPARSE_FILL_MODE_FULL; break;
	case Matrix_Shape::UPPER: descr.mode = SPARSE_FILL_MODE_UPPER; break;
	case Matrix_Shape::LOWER: descr.mode = SPARSE_FILL_MODE_LOWER; break;
	}
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_IJV<double> &a, const double &beta, const Vector_Dense<double> &x)
{
	matrix_descr descr;
	switch (a.type) {
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

	switch (a.shape) {
	case Matrix_Shape::FULL : descr.mode = SPARSE_FILL_MODE_FULL; break;
	case Matrix_Shape::UPPER: descr.mode = SPARSE_FILL_MODE_UPPER; break;
	case Matrix_Shape::LOWER: descr.mode = SPARSE_FILL_MODE_LOWER; break;
	}
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._spblas->inspector, descr, x.vals, beta, y.vals));
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_IJV<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	matrix_descr descr;
	switch (a.type) {
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

	switch (a.shape) {
	case Matrix_Shape::FULL : descr.mode = SPARSE_FILL_MODE_FULL; break;
	case Matrix_Shape::UPPER: descr.mode = SPARSE_FILL_MODE_UPPER; break;
	case Matrix_Shape::LOWER: descr.mode = SPARSE_FILL_MODE_LOWER; break;
	}
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._spblas->inspector, descr, x.vals, beta, y.vals));
}

}
}

#endif
#endif
