
#include "math/math.h"
#include "esinfo/eslog.h"
#include <cstdio>

#ifdef HAVE_MKL
#include "w.mkl.h"
#include "mkl_spblas.h"
#endif

using namespace espreso;

#ifdef HAVE_MKL
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
#endif

#ifdef HAVE_MKL
void MATH::CSRHandlerData::info(esint &nrows, esint &ncols, esint &nnz)
{
	sparse_index_base_t indexing;
	esint *rstart, *rend, *cids;
	double *vals;
	checkStatus(mkl_sparse_d_export_csr(inspector, &indexing, &nrows, &ncols, &rstart, &rend, &cids, &vals));
	nnz = rend[nrows - 1] - rstart[0];
}

void MATH::CSRHandlerData::info(esint &nrows, esint &ncols, esint &nnz, esint* &rows, esint* &cols, double* &vals)
{
	sparse_index_base_t indexing;
	esint *rend;
	checkStatus(mkl_sparse_d_export_csr(inspector, &indexing, &nrows, &ncols, &rows, &rend, &cols, &vals));
	nnz = rows[nrows] - rows[0];
}
#endif

void MATH::CSRHandler::sizes(esint &nrows, esint &ncols, esint &nnz) const
{
#ifdef HAVE_MKL
	inner->info(nrows, ncols, nnz);
#endif
}

MATH::CSRHandler::CSRHandler(esint nrows, esint ncols, esint nnz, esint *rows, esint *cols, double *vals)
{
#ifdef HAVE_MKL
	inner = new CSRHandlerData();
	if (nnz) {
		checkStatus(mkl_sparse_d_create_csr(&inner->inspector, SPARSE_INDEX_BASE_ONE, nrows, ncols, rows, rows + 1, cols, vals));
	}
#endif
}

void MATH::CSRHandler::data(esint* &rows, esint* &cols, double* &vals) const
{
#ifdef HAVE_MKL
	esint nrows, ncols, nnz;
	inner->info(nrows, ncols, nnz, rows, cols, vals);
#endif
}

MATH::CSRHandler::~CSRHandler()
{
#ifdef HAVE_MKL
	if (inner->called) {
		CSRMatClearFactors(this);
	}
	if (inner->inspector) {
		checkStatus(mkl_sparse_destroy(inner->inspector));
	}
	delete inner;
#endif
}

MATH::IJVHandler::IJVHandler(esint nrows, esint ncols, esint nnz, esint *rows, esint *cols, double *vals)
{
#ifdef HAVE_MKL
	inner = new IJVHandlerData();
	if (nnz) { // empty handler is not allowed
		checkStatus(mkl_sparse_d_create_coo(&inner->inspector, SPARSE_INDEX_BASE_ONE, nrows, ncols, nnz, rows, cols, vals));
	}
#endif
}

MATH::IJVHandler::~IJVHandler()
{
#ifdef HAVE_MKL
	if (inner->inspector) {
		checkStatus(mkl_sparse_destroy(inner->inspector));
	}
	delete inner;
#endif
}

void MATH::upCSRMatVecProduct(CSRHandler *A, double *vals, double *result)
{
#ifdef HAVE_MKL
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A->inner->inspector, descr, vals, 0, result));
#endif
}

void MATH::upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, float *mVals, float *vVals, float *result)
{
#ifdef HAVE_MKL
	sparse_matrix_t inspector;
	checkStatus(mkl_sparse_s_create_csr(&inspector, SPARSE_INDEX_BASE_ONE, rows, cols, mRows, mRows + 1, mCols, mVals));
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, inspector, descr, vVals, 0, result));
	checkStatus(mkl_sparse_destroy(inspector));
#endif
}

void MATH::upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, double *mVals, double *vVals, double *result)
{
#ifdef HAVE_MKL
	sparse_matrix_t inspector;
	checkStatus(mkl_sparse_d_create_csr(&inspector, SPARSE_INDEX_BASE_ONE, rows, cols, mRows, mRows + 1, mCols, mVals));
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, inspector, descr, vVals, 0, result));
	checkStatus(mkl_sparse_destroy(inspector));
#endif
}

void MATH::CSRMatVecProduct(CSRHandler *A, double *vals, double *result)
{
#ifdef HAVE_MKL
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	descr.mode = SPARSE_FILL_MODE_FULL;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A->inner->inspector, descr, vals, 0, result));
#endif
}

void MATH::CSRMatCSRMatProduct(CSRHandler *C, CSRHandler *A, CSRHandler *B, bool transposeA)
{
#ifdef HAVE_MKL
	if (C->inner->inspector) {
		checkStatus(mkl_sparse_destroy(C->inner->inspector));
	}
	if (transposeA) {
		checkStatus(mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, A->inner->inspector, B->inner->inspector, &C->inner->inspector));
	} else {
		checkStatus(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A->inner->inspector, B->inner->inspector, &C->inner->inspector));
	}
	checkStatus(mkl_sparse_order(C->inner->inspector));
#endif
}

void MATH::CSRTranspose(CSRHandler *A, CSRHandler *At)
{
#ifdef HAVE_MKL
	eslog::internalFailure("cannot use 'mkl_sparse_convert_csr'. It has weird behaviour.\n");
	checkStatus(mkl_sparse_convert_csr(A->inner->inspector, SPARSE_OPERATION_TRANSPOSE, &At->inner->inspector));
#endif
}

