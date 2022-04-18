
#include "math/math.h"
#include "math/math.h"
#include "esinfo/eslog.h"
#include <cstdio>

#ifdef HAVE_MKL
#include "w.mkl.h"
#include "mkl_spblas.h"
#include "mkl.h"
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

template<typename T>
int _pardisoType(const Matrix_CSR<T> &x)
{
	switch (x.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:    return  2;
	case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:           return -2;
	case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:         return  1;
	case Matrix_Type::REAL_NONSYMMETRIC:                   return 11;
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE: return  4;
	case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:        return -4;
	case Matrix_Type::COMPLEX_SYMMETRIC:                   return  6;
	case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:      return  3;
	case Matrix_Type::COMPLEX_NONSYMMETRIC:                return 13;
	}
	return 0;
}

template<typename T>
bool _callPardiso(esint phase, const Matrix_CSR<T> &m, esint nrhs, T *rhs, T *solution)
{
	m._external->phase = phase;
	pardiso(
			m._external->pt, &m._external->maxfct, &m._external->mnum,
			&m._external->mtype,
			&m._external->phase,
			&m.nrows, m.vals, m.rows, m.cols,
			m._external->perm, &nrhs, m._external->iparm, &m._external->msglvl,
			rhs, solution,
			&m._external->error);

	switch (m._external->error) {
	case   0: break;
	case  -1: eslog::error("MKL PARDISO: input inconsistent.\n"); break;
	case  -2: eslog::error("MKL PARDISO: not enough memory.\n"); break;
	case  -3: eslog::error("MKL PARDISO: reordering problem.\n"); break;
	case  -4: eslog::error("MKL PARDISO: zero pivot, numerical factorization or iterative refinement problem.\n"); break;
	case  -5: eslog::error("MKL PARDISO: unclassified (internal) error.\n"); break;
	case  -6: eslog::error("MKL PARDISO: reordering failed.\n"); break;
	case  -7: eslog::error("MKL PARDISO: diagonal matrix is singular.\n"); break;
	case  -8: eslog::error("MKL PARDISO: 32-bit integer overflow problem.\n"); break;
	case  -9: eslog::error("MKL PARDISO: not enough memory for OOC.\n"); break;
	case -10: eslog::error("MKL PARDISO: error opening OOC files.\n"); break;
	case -11: eslog::error("MKL PARDISO: read/write error with OOC files.\n"); break;
	case -12: eslog::error("MKL PARDISO: (pardiso_64 only) pardiso_64 called from 32-bit library.\n"); break;
	case -13: eslog::error("MKL PARDISO: interrupted by the (user-defined) mkl_progress function.\n"); break;
	case -15: eslog::error("MKL PARDISO: internal error which can appear for iparm[23]=10 and iparm[12]=1. Try switch matching off (set iparm[12]=0 and rerun.\n"); break;
	}
	return m._external->error == 0;
}
#endif

namespace espreso {
namespace math {

template <>
void commit(Matrix_Dense<double> &x)
{
#ifdef HAVE_MKL
	// dummy
#endif
}

template <>
void commit(Matrix_Dense<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	// dummy
#endif
}

template <>
void commit(Matrix_CSR<double> &x)
{
#ifdef HAVE_MKL
	x._external = new Matrix_CSR_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_d_create_csr(&x._external->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.rows, x.rows + 1, x.cols, x.vals));
	}
	x._external->mtype = _pardisoType(x);
	pardisoinit(x._external->pt, &x._external->mtype, x._external->iparm);
	x._external->iparm[0] = 1;			/* No solver default */
	x._external->iparm[1] = 2;			/* Fill-in reordering from METIS */
	x._external->iparm[2] = 1; 			//by default the solver runs with single thread
	x._external->iparm[9] = 13;			/* Perturb the pivot elements with 1E-13 */
	x._external->iparm[10] = 1;			/* Use nonsymmetric permutation and scaling MPS */
#endif
}

template <>
void commit(Matrix_CSR<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	x._external = new Matrix_CSR_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_z_create_csr(&x._external->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.rows, x.rows + 1, x.cols, x.vals));
	}
	x._external->mtype = _pardisoType(x);
	pardisoinit(x._external->pt, &x._external->mtype, x._external->iparm);
	x._external->iparm[0] = 1;			/* No solver default */
	x._external->iparm[1] = 2;			/* Fill-in reordering from METIS */
	x._external->iparm[2] = 1; 			//by default the solver runs with single thread
	x._external->iparm[9] = 13;			/* Perturb the pivot elements with 1E-13 */
	x._external->iparm[10] = 1;			/* Use nonsymmetric permutation and scaling MPS */
#endif
}

template <>
void commit(Matrix_IJV<double> &x)
{
#ifdef HAVE_MKL
	x._external = new Matrix_IJV_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_d_create_coo(&x._external->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.nnz, x.rows, x.cols, x.vals));
	}
#endif
}

template <>
void commit(Matrix_IJV<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	x._external = new Matrix_IJV_External_Representation();
	if (x.nnz) {
		checkStatus(mkl_sparse_z_create_coo(&x._external->inspector, SPARSE_INDEX_BASE_ONE, x.nrows, x.ncols, x.nnz, x.rows, x.cols, x.vals));
	}
#endif
}

template <>
void free(Matrix_Dense<double> &x)
{
#ifdef HAVE_MKL
#endif
}

template <>
void free(Matrix_Dense<std::complex<double> > &x)
{
#ifdef HAVE_MKL
#endif
}

template <>
void free(Matrix_CSR<double> &x)
{
#ifdef HAVE_MKL
	if (x._external) {
		checkStatus(mkl_sparse_destroy(x._external->inspector));
		x._external = nullptr;
		_callPardiso<double>(-1, x, 0, nullptr, nullptr);
	}
#endif
}

template <>
void free(Matrix_CSR<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	if (x._external) {
		checkStatus(mkl_sparse_destroy(x._external->inspector));
		x._external = nullptr;
		_callPardiso<std::complex<double> >(-1, x, 0, nullptr, nullptr);
	}
#endif
}

template <>
void free(Matrix_IJV<double> &x)
{
#ifdef HAVE_MKL
	if (x._external) {
		checkStatus(mkl_sparse_destroy(x._external->inspector));
		x._external = nullptr;
	}
#endif
}

template <>
void free(Matrix_IJV<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	if (x._external) {
		checkStatus(mkl_sparse_destroy(x._external->inspector));
		x._external = nullptr;
	}
#endif
}

template <>
void symbolicFactorization(const Matrix_CSR<double> &x)
{
#ifdef HAVE_MKL
	_callPardiso<double>(11, x, 0, nullptr, nullptr);
#endif
}

template <>
void symbolicFactorization(const Matrix_CSR<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	_callPardiso<std::complex<double> >(11, x, 0, nullptr, nullptr);
#endif
}

template <>
void numericalFactorization(const Matrix_CSR<double> &x)
{
#ifdef HAVE_MKL
	_callPardiso<double>(22, x, 0, nullptr, nullptr);
#endif
}

template <>
void numericalFactorization(const Matrix_CSR<std::complex<double> > &x)
{
#ifdef HAVE_MKL
	_callPardiso<std::complex<double> >(22, x, 0, nullptr, nullptr);
#endif
}

template <>
void solve(const Matrix_CSR<double> &x, Vector_Dense<double> &rhs, Vector_Dense<double> &solution)
{
#ifdef HAVE_MKL
	_callPardiso<double>(33, x, 1, rhs.vals, solution.vals);
#endif
}

template <>
void solve(const Matrix_CSR<double> &x, Matrix_Dense<double> &rhs, Matrix_Dense<double> &solution)
{
#ifdef HAVE_MKL
	_callPardiso<double>(33, x, rhs.nrows, rhs.vals, solution.vals);
#endif
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &x, Vector_Dense<std::complex<double> > &rhs, Vector_Dense<std::complex<double> > &solution)
{
#ifdef HAVE_MKL
	_callPardiso<std::complex<double> >(33, x, 1, rhs.vals, solution.vals);
#endif
}

template <>
void solve(const Matrix_CSR<std::complex<double> > &x, Matrix_Dense<std::complex<double> > &rhs, Matrix_Dense<std::complex<double> > &solution)
{
#ifdef HAVE_MKL
	_callPardiso<std::complex<double> >(33, x, rhs.nrows, rhs.vals, solution.vals);
#endif
}

template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_CSR<double> &a, const double &beta, const Vector_Dense<double> &x)
{
#ifdef HAVE_MKL
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
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._external->inspector, descr, x.vals, beta, y.vals));
#endif
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_CSR<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
#ifdef HAVE_MKL
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
	checkStatus(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._external->inspector, descr, x.vals, beta, y.vals));
#endif
}

template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_IJV<double> &a, const double &beta, const Vector_Dense<double> &x)
{
#ifdef HAVE_MKL
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
	checkStatus(mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._external->inspector, descr, x.vals, beta, y.vals));
#endif
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_IJV<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
#ifdef HAVE_MKL
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
	checkStatus(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, a._external->inspector, descr, x.vals, beta, y.vals));
#endif
}

}
}
