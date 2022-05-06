
#ifndef SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_CHOLMOD_H_
#define SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_CHOLMOD_H_

#include "cholmod.h"
#include "math/wrappers/math.spblas.h"
#include "math/wrappers/math.solver.h"

#include <cstring>

namespace espreso {

template <typename T> inline void _start(cholmod_common &common);
template <> inline void _start<int>(cholmod_common &common) { cholmod_start(&common); }
template <> inline void _start<long>(cholmod_common &common) { cholmod_l_start(&common); }

template <typename T> inline void _analyze(cholmod_factor* &L, cholmod_sparse* A, cholmod_common &common);
template <> inline void _analyze<int>(cholmod_factor* &L, cholmod_sparse* A, cholmod_common &common) { A->itype = CHOLMOD_INT; L = cholmod_analyze(A, &common); }
template <> inline void _analyze<long>(cholmod_factor* &L, cholmod_sparse* A, cholmod_common &common) { A->itype = CHOLMOD_LONG; L = cholmod_l_analyze(A, &common); }

template <typename T> inline void _factorize(cholmod_factor *L, cholmod_sparse *A, cholmod_common &common);
template <> inline void _factorize<int>(cholmod_factor *L, cholmod_sparse *A, cholmod_common &common) { cholmod_factorize(A, L, &common); }
template <> inline void _factorize<long>(cholmod_factor *L, cholmod_sparse *A, cholmod_common &common) { cholmod_l_factorize(A, L, &common); }

template <typename T> inline void _solve(cholmod_dense* &x, cholmod_factor *L, cholmod_dense *b, cholmod_common &common);
template <> inline void _solve<int>(cholmod_dense* &x, cholmod_factor *L, cholmod_dense *b, cholmod_common &common) { x = cholmod_solve(CHOLMOD_A, L, b, &common); }
template <> inline void _solve<long>(cholmod_dense* &x, cholmod_factor *L, cholmod_dense *b, cholmod_common &common) { x = cholmod_l_solve(CHOLMOD_A, L, b, &common); }

template <typename T> inline void _apply(cholmod_dense* Y, cholmod_sparse* A, cholmod_dense *X, double alpha[2], double beta[2], cholmod_common &common);
template <> inline void _apply<int>(cholmod_dense* Y, cholmod_sparse* A, cholmod_dense *X, double alpha[2], double beta[2], cholmod_common &common) { cholmod_sdmult(A, 0, alpha, beta, X, Y, &common); }
template <> inline void _apply<long>(cholmod_dense* Y, cholmod_sparse* A, cholmod_dense *X, double alpha[2], double beta[2], cholmod_common &common) { cholmod_sdmult(A, 0, alpha, beta, X, Y, &common); }

template <typename T> inline void _finish(cholmod_common &common);
template <> inline void _finish<int>(cholmod_common &common) { cholmod_finish(&common); }
template <> inline void _finish<long>(cholmod_common &common) { cholmod_l_finish(&common); }

template <typename T, typename O> inline void _free(O* object, cholmod_common &common);
template <> inline void _free<int, cholmod_sparse>(cholmod_sparse* object, cholmod_common &common) { cholmod_free_sparse(&object, &common); }
template <> inline void _free<long, cholmod_sparse>(cholmod_sparse* object, cholmod_common &common) { cholmod_l_free_sparse(&object, &common); }
template <> inline void _free<int, cholmod_dense>(cholmod_dense* object, cholmod_common &common) { cholmod_free_dense(&object, &common); }
template <> inline void _free<long, cholmod_dense>(cholmod_dense* object, cholmod_common &common) { cholmod_l_free_dense(&object, &common); }
template <> inline void _free<int, cholmod_factor>(cholmod_factor* object, cholmod_common &common) { cholmod_free_factor(&object, &common); }
template <> inline void _free<long, cholmod_factor>(cholmod_factor* object, cholmod_common &common) { cholmod_l_free_factor(&object, &common); }


template <typename T>
static inline void set(cholmod_sparse *A, const Matrix_CSR<T> &M)
{
	A->nrow = M.nrows;
	A->ncol = M.ncols;
	A->nzmax = M.nnz;

	A->p = M.rows;
	A->i = M.cols;
	A->x = M.vals;

	A->stype = -1; // CSR -> CSC
	A->xtype = CHOLMOD_PATTERN;
	A->dtype = CHOLMOD_DOUBLE;
	A->sorted = 1;
	A->packed = 1;
}

static inline void update(cholmod_sparse *A, const Matrix_CSR<double> &M)
{
	A->x = M.vals;
	A->xtype = CHOLMOD_REAL;
}

static inline void update(cholmod_sparse *A, const Matrix_CSR<std::complex<double> > &M)
{
	A->x = M.vals;
	A->xtype = CHOLMOD_COMPLEX;
}

static inline void update(cholmod_dense *A, const Matrix_Dense<double> &M)
{
	A->nrow = M.ncols;
	A->ncol = M.nrows;
	A->nzmax = M.ncols * M.nrows;
	A->d = M.ncols;
	A->x = M.vals;
	A->xtype = CHOLMOD_REAL;
	A->dtype = CHOLMOD_DOUBLE;
}

static inline void update(cholmod_dense *A, const Matrix_Dense<std::complex<double> > &M)
{
	A->nrow = M.nrows;
	A->ncol = M.ncols;
	A->nzmax = M.ncols * M.nrows;
	A->d = M.nrows;
	A->x = M.vals;
	A->xtype = CHOLMOD_COMPLEX;
	A->dtype = CHOLMOD_DOUBLE;
}

static inline void update(cholmod_dense *A, const Vector_Dense<double> &v)
{
	A->nrow = v.size;
	A->ncol = 1;
	A->nzmax = v.size;
	A->d = v.size;
	A->x = v.vals;
	A->xtype = CHOLMOD_REAL;
	A->dtype = CHOLMOD_DOUBLE;
}

static inline void update(cholmod_dense *A, const Vector_Dense<std::complex<double> > &v)
{
	A->nrow = v.size;
	A->ncol = 1;
	A->nzmax = v.size;
	A->d = v.size;
	A->x = v.vals;
	A->xtype = CHOLMOD_COMPLEX;
	A->dtype = CHOLMOD_DOUBLE;
}

static inline void extract(cholmod_dense *A, cholmod_common &common, Matrix_Dense<double> &M)
{
	memcpy(M.vals, A->x, sizeof(double) * M.nrows * M.ncols);
	_free<esint>(A, common);
}

static inline void extract(cholmod_dense *A, cholmod_common &common, Matrix_Dense<std::complex<double> > &M)
{
	memcpy(M.vals, A->x, sizeof(std::complex<double>) * M.nrows * M.ncols);
	_free<esint>(A, common);
}

static inline void extract(cholmod_dense *A, cholmod_common &common, Vector_Dense<double> &v)
{
	memcpy(v.vals, A->x, sizeof(double) * v.size);
	_free<esint>(A, common);
}

static inline void extract(cholmod_dense *A, cholmod_common &common, Vector_Dense<std::complex<double> > &v)
{
	memcpy(v.vals, A->x, sizeof(std::complex<double>) * v.size);
	_free<esint>(A, common);
}

}

#endif /* SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_CHOLMOD_H_ */
