
#include "w.mkl.pdss.h"
#include "esinfo/mpiinfo.h"
#include "math2/primitives/matrix_csr.h"

#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"

#include <vector>
#include <complex>

#ifdef HAVE_MKLPDSS
#include "mkl_cluster_sparse_solver.h"

namespace espreso {

template<typename T>
struct MKLPDSSDataHolder {
	void *pt[64];

	esint maxfct;
	esint mnum;

	esint mtype;
	esint phase;
	esint perm;
	esint n;
	esint nrhs;
	esint iparm[64];
	esint msglvl;
	int   comm;
	esint error;

	Matrix_CSR<T> A;
	Vector_Dense<T> b, x;
//
//	esint *rowPtrs;
//	esint *colIndices;
//	double *values;
//	double *rhsValues;
//	double *solution;
};

}

#endif

namespace espreso {

void _check()
{
#ifndef HAVE_MKLPDSS
	eslog::globalerror("ESPRESO run-time error: cannot call MKL PDSS solver (the library with the solver is not linked).\n");
#endif
}

template<typename T>
bool _call(MKLPDSS<T> &mklpdss, esint phase)
{
#ifdef HAVE_MKLPDSS
	mklpdss.external->phase = phase;
	cluster_sparse_solver(
			mklpdss.external->pt, &mklpdss.external->maxfct, &mklpdss.external->mnum,
			&mklpdss.external->mtype,
			&mklpdss.external->phase,
			&mklpdss.external->n, mklpdss.external->A.vals, mklpdss.external->A.rows, mklpdss.external->A.cols,
			&mklpdss.external->perm, &mklpdss.external->nrhs, mklpdss.external->iparm, &mklpdss.external->msglvl,
			mklpdss.external->b.vals, mklpdss.external->x.vals,
			&mklpdss.external->comm, &mklpdss.external->error);

	switch (mklpdss.external->error) {
	case   0: break;
	case  -1: eslog::error("MKL PDSS: input inconsistent.\n"); break;
	case  -2: eslog::error("MKL PDSS: not enough memory.\n"); break;
	case  -3: eslog::error("MKL PDSS: reordering problem.\n"); break;
	case  -4: eslog::error("MKL PDSS: zero pivot, numerical factorization or iterative refinement problem.\n"); break;
	case  -5: eslog::error("MKL PDSS: unclassified (internal) error.\n"); break;
	case  -6: eslog::error("MKL PDSS: reordering failed.\n"); break;
	case  -7: eslog::error("MKL PDSS: diagonal matrix is singular.\n"); break;
	case  -8: eslog::error("MKL PDSS: 32-bit integer overflow problem.\n"); break;
	case  -9: eslog::error("MKL PDSS: not enough memory for OOC.\n"); break;
	case -10: eslog::error("MKL PDSS: error opening OOC files.\n"); break;
	case -11: eslog::error("MKL PDSS: read/write error with OOC files.\n"); break;
	}
	return mklpdss.external->error == 0;
#endif
}

template<typename T>
bool _set(MKLPDSS<T> &mklpdss, const Matrix_Distributed<Matrix_CSR, T> &A)
{
#ifdef HAVE_MKLPDSS
	double start = eslog::time();

	mklpdss.external = new MKLPDSSDataHolder<T>();
	mklpdss.external->n = A.distribution.totalSize;

	mklpdss.external->maxfct = 1; // dummy
	mklpdss.external->mnum = 1; // dummy

	mklpdss.external->nrhs = 1;

	std::fill(mklpdss.external->iparm, mklpdss.external->iparm + 64, 0);
	mklpdss.external->iparm[0] = 1; // Use filled values.
	// Fill-in reducing ordering for the input matrix.
	mklpdss.external->iparm[1] = info::mpi::size > 1 ? 10 : 3; // MPI or parallel
	// Matrix input format.
	mklpdss.external->iparm[39] = 2; // distributed A, x, rhs
	mklpdss.external->iparm[40] = A.distribution.begin + 1;
	mklpdss.external->iparm[41] = A.distribution.end;

	mklpdss.external->msglvl = 0;
	mklpdss.external->comm = MPI_Comm_c2f(info::mpi::comm);

	switch (A.cluster.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		eslog::solver("     - | MATRIX TYPE ::                               REAL SYMMETRIC POSITIVE DEFINITE | -\n");
		mklpdss.external->mtype = 2; break;
	case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:
		eslog::solver("     - | MATRIX TYPE ::                                      REAL SYMMETRIC INDEFINITE | -\n");
		mklpdss.external->mtype = -2; break;
	case Matrix_Type::REAL_UNSYMMETRIC:
		eslog::solver("     - | MATRIX TYPE ::                      REAL NONSYMMETRIC (STRUCTURALY SYMMETRIC) | -\n");
		mklpdss.external->mtype = 1; break;
	}

	// pick only upper triangle (since composer does not set correct dirichlet in symmetric matrices)
	if (A.cluster.type == Matrix_Type::REAL_SYMMETRIC_INDEFINITE || A.cluster.type == Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		esint nhalo = A.distribution.halo.size();
		for (esint i = nhalo; i < A.cluster.nrows; i++) {
			for (esint c = A.cluster.rows[i] - _Matrix_CSR_Pattern::Indexing; c < A.cluster.rows[i + 1] - _Matrix_CSR_Pattern::Indexing; ++c) {
				if (A.distribution.begin + i - nhalo <= A.cluster.cols[c] - _Matrix_CSR_Pattern::Indexing) {
					++mklpdss.external->A.nnz;
				}
			}
		}
		mklpdss.external->A.resize(A.cluster.nrows - nhalo, A.cluster.ncols, mklpdss.external->A.nnz);
		mklpdss.external->A.rows[0] = _Matrix_CSR_Pattern::Indexing;
		for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
			for (esint c = A.cluster.rows[i] - _Matrix_CSR_Pattern::Indexing; c < A.cluster.rows[i + 1] - _Matrix_CSR_Pattern::Indexing; ++c) {
				if (A.distribution.begin + i - nhalo <= A.cluster.cols[c] - _Matrix_CSR_Pattern::Indexing) {
					mklpdss.external->A.cols[offset++] = A.cluster.cols[c];
				}
			}
			mklpdss.external->A.rows[i - nhalo + 1] = offset + _Matrix_CSR_Pattern::Indexing;
		}
	} else {
		esint nhalo = A.distribution.halo.size();
		for (esint i = nhalo; i < A.cluster.nrows; i++) {
			for (esint c = A.cluster.rows[i] - _Matrix_CSR_Pattern::Indexing; c < A.cluster.rows[i + 1] - _Matrix_CSR_Pattern::Indexing; ++c) {
				++mklpdss.external->A.nnz;
			}
		}
		mklpdss.external->A.resize(A.cluster.nrows - nhalo, A.cluster.ncols, mklpdss.external->A.nnz);
		mklpdss.external->A.rows[0] = _Matrix_CSR_Pattern::Indexing;
		for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
			for (esint c = A.cluster.rows[i] - _Matrix_CSR_Pattern::Indexing; c < A.cluster.rows[i + 1] - _Matrix_CSR_Pattern::Indexing; ++c) {
				mklpdss.external->A.cols[offset++] = A.cluster.cols[c];
			}
			mklpdss.external->A.rows[i - nhalo + 1] = offset + _Matrix_CSR_Pattern::Indexing;
		}
	}
	bool status = _call(mklpdss, 11);
	eslog::solver("     - | SYMBOLIC FACTORIZATION                                             %8.3f s | -\n", eslog::time() - start);
	return status;
#endif
}

template<typename T>
bool _update(MKLPDSS<T> &mklpdss, const Matrix_Distributed<Matrix_CSR, T> &A)
{
#ifdef HAVE_MKLPDSS
	double start = eslog::time();
	if (A.cluster.type == Matrix_Type::REAL_SYMMETRIC_INDEFINITE || A.cluster.type == Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		esint nhalo = A.distribution.halo.size();
		for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
			for (esint c = A.cluster.rows[i] - _Matrix_CSR_Pattern::Indexing; c < A.cluster.rows[i + 1] - _Matrix_CSR_Pattern::Indexing; ++c) {
				if (A.distribution.begin + i - nhalo <= A.cluster.cols[c] - _Matrix_CSR_Pattern::Indexing) {
					mklpdss.external->A.vals[offset++] = A.cluster.vals[c];
				}
			}
		}
	} else {
		esint nhalo = A.distribution.halo.size();
		for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
			for (esint c = A.cluster.rows[i] - _Matrix_CSR_Pattern::Indexing; c < A.cluster.rows[i + 1] - _Matrix_CSR_Pattern::Indexing; ++c) {
				mklpdss.external->A.vals[offset++] = A.cluster.vals[c];
			}
		}
	}
	bool status = _call(mklpdss, 22);
	eslog::solver("     - | NUMERICAL FACTORIZATION                                            %8.3f s | -\n", eslog::time() - start);
	return status;
#endif
}

template<typename T>
bool _solve(MKLPDSS<T> &mklpdss, const Vector_Distributed<Vector_Dense, T> &b, Vector_Distributed<Vector_Dense, T> &x)
{
#ifdef HAVE_MKLPDSS
	mklpdss.external->b.vals = b.cluster.vals + b.distribution.halo.size();
	mklpdss.external->x.vals = x.cluster.vals + x.distribution.halo.size();
	double start = eslog::time();

	bool status = _call(mklpdss, 33); // solve at once
//	_data.x.scatterToUpper();
	eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
	eslog::solver("     - --------------------------------------------------------------------------------- -\n");
	return status;
#endif
}

template<typename T>
void _clear(MKLPDSS<T> &mklpdss)
{
#ifdef HAVE_MKLPDSS
	delete mklpdss.external;
#endif
}

template<> bool MKLPDSS<double>::set(const Matrix_Distributed<Matrix_CSR, double> &A) { return _set(*this, A); }
template<> bool MKLPDSS<std::complex<double> >::set(const Matrix_Distributed<Matrix_CSR, std::complex<double> > &A) { return _set(*this, A); }

template<> bool MKLPDSS<double>::update(const Matrix_Distributed<Matrix_CSR, double> &A) { return _update(*this, A); }
template<> bool MKLPDSS<std::complex<double> >::update(const Matrix_Distributed<Matrix_CSR, std::complex<double> > &A) { return _update(*this, A); }

template<> bool MKLPDSS<double>::solve(const Vector_Distributed<Vector_Dense, double> &b, Vector_Distributed<Vector_Dense, double> &x) { return _solve(*this, b, x); }
template<> bool MKLPDSS<std::complex<double> >::solve(const Vector_Distributed<Vector_Dense, std::complex<double> > &b, Vector_Distributed<Vector_Dense, std::complex<double>> &x) { return _solve(*this, b, x); }

template<> void MKLPDSS<double>::check() { _check(); }
template<> void MKLPDSS<std::complex<double> >::check() { _check(); }

template<> bool MKLPDSS<double>::call(int phase) { return _call(*this, phase); }
template<> bool MKLPDSS<std::complex<double> >::call(int phase) { return _call(*this, phase); }

template<> void MKLPDSS<double>::clear() { _clear(*this); }
template<> void MKLPDSS<std::complex<double> >::clear() { _clear(*this); }

}
