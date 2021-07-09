
#include "w.mkl.pdss.systemsolver.h"
#include "physics/system/mklpdsssystem.h"
#include "esinfo/mpiinfo.h"

#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"

#include <vector>
#include <complex>

#ifdef HAVE_MKLPDSS
#include "mkl_cluster_sparse_solver.h"

namespace espreso {
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

	std::vector<esint> rowData, colData;
	std::vector<double> valData;

	esint *rowPtrs;
	esint *colIndices;
	double *values;
	double *rhsValues;
	double *solution;
};
}
#endif

using namespace espreso;

MKLPDSSSystemSolver::MKLPDSSSystemSolver(MKLPDSSConfiguration &configuration, MKLPDSSSolverData &data)
: configuration(configuration), _roffset(0), _nrows(0), _precision(1e-15), _data(data), _inner(NULL)
{
#ifndef HAVE_MKLPDSS
	eslog::globalerror("ESPRESO run-time error: cannot call MKL PDSS solver (the library with the solver is not linked).\n");
#endif
}

void MKLPDSSSystemSolver::init()
{
#ifdef HAVE_MKLPDSS
	_nrows = _data.K.nrows - _data.K.nhalo;
	_roffset = _nrows;

	_inner = new MKLPDSSDataHolder();
	_inner->n = Communication::exscan(_roffset);

	_inner->maxfct = 1; // dummy
	_inner->mnum = 1; // dummy

	_inner->nrhs = 1;

	std::fill(_inner->iparm, _inner->iparm + 64, 0);
	_inner->iparm[0] = 1; // Use filled values.
	// Fill-in reducing ordering for the input matrix.
	_inner->iparm[1] = info::mpi::size > 1 ? 10 : 3; // MPI or parallel
	// Matrix input format.
	_inner->iparm[39] = 2; // distributed A, x, rhs
	_inner->iparm[40] = _roffset + 1;
	_inner->iparm[41] = _roffset + _nrows;

	_inner->msglvl = 0;
	_inner->comm = MPI_Comm_c2f(info::mpi::comm);
#endif
}

bool MKLPDSSSystemSolver::update()
{
#ifdef HAVE_MKLPDSS
	eslog::solver("     - ---- LINEAR SOLVER -------------------------------------------------------------- -\n");
	eslog::solver("     - | SOLVER ::     MKL                       TYPE :: PARALLEL DIRECT SPARSE SOLVER | -\n");
	double start = eslog::time();

	switch (_data.K.type) {
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		eslog::solver("     - | MATRIX TYPE ::                               REAL SYMMETRIC POSITIVE DEFINITE | -\n");
		_inner->mtype = 2; break;
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
		eslog::solver("     - | MATRIX TYPE ::                                      REAL SYMMETRIC INDEFINITE | -\n");
		_inner->mtype = -2; break;
	case MatrixType::REAL_UNSYMMETRIC:
		eslog::solver("     - | MATRIX TYPE ::                      REAL NONSYMMETRIC (STRUCTURALY SYMMETRIC) | -\n");
		_inner->mtype = 1; break;
	}

	const bool fillStructure = !_inner->rowData.size();

	// pick only upper triangle (since composer does not set correct dirichlet in symmetric matrices)
	if (_data.K.type == MatrixType::REAL_SYMMETRIC_INDEFINITE || _data.K.type == MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		if (fillStructure) {
			_inner->rowData.clear();
			_inner->colData.clear();
			_inner->rowData.reserve(_nrows + 1);
			_inner->rowData.push_back(1);
		}
		_inner->valData.clear();
		for (esint i = _data.K.nhalo; i < _data.K.nrows; i++) {
			for (esint c = _data.K.rows[i] - 1; c < _data.K.rows[i + 1] - 1; c++) {
				if (_roffset + i - _data.K.nhalo <= _data.K.cols[c] - 1) {
					if (fillStructure) {
						_inner->colData.push_back(_data.K.cols[c]);
					}
					_inner->valData.push_back(_data.K.vals[c]);
				}
			}
			if (fillStructure) {
				_inner->rowData.push_back(_inner->colData.size() + 1);
			}
		}
		_inner->colIndices = _inner->colData.data();
		_inner->values = _inner->valData.data();
	} else {
		if (fillStructure) {
			_inner->rowData.clear();
			_inner->rowData.reserve(_nrows + 1);

			// row data have to be always renumbered
			for (esint i = _data.K.nhalo; i <= _data.K.nrows; i++) {
				_inner->rowData.push_back(_data.K.rows[i] - _data.K.rows[_data.K.nhalo] + 1);
			}
			_inner->colIndices = _data.K.cols + _data.K.rows[_data.K.nhalo] - 1;
		}
		_inner->values = _data.K.vals + _data.K.rows[_data.K.nhalo] - 1;
	}

	_inner->rowPtrs = _inner->rowData.data();

	if (fillStructure) {
		if (!call(11)) return false;
		eslog::solver("     - | SYMBOLIC FACTORIZATION                                             %8.3f s | -\n", eslog::time() - start);
		start = eslog::time();
	}
	if (!call(22)) return false;
	eslog::solver("     - | NUMERICAL FACTORIZATION                                            %8.3f s | -\n", eslog::time() - start);

	_inner->rhsValues = _data.f[0].vals + _data.f[0].nhalo;
#endif
	return true;
}

bool MKLPDSSSystemSolver::solve()
{
#ifdef HAVE_MKLPDSS
	_inner->solution = _data.x[0].vals + _data.x[0].nhalo;
	double start = eslog::time();
	if (!call(33)) return false; // solve at once
	_data.x.scatterToUpper();
	eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
	eslog::solver("     - --------------------------------------------------------------------------------- -\n");

	return true;
#endif
}

MKLPDSSSystemSolver::~MKLPDSSSystemSolver()
{
#ifdef HAVE_MKLPDSS
	call(-1);
	delete _inner;
#endif
}

bool MKLPDSSSystemSolver::call(esint phase)
{
#ifdef HAVE_MKLPDSS
	_inner->phase = phase;
	cluster_sparse_solver(
			_inner->pt, &_inner->maxfct, &_inner->mnum,
			&_inner->mtype,
			&_inner->phase,
			&_inner->n, _inner->values, _inner->rowPtrs, _inner->colIndices,
			&_inner->perm, &_inner->nrhs, _inner->iparm, &_inner->msglvl,
			_inner->rhsValues, _inner->solution,
			&_inner->comm, &_inner->error);

	bool ret = true;
	switch (_inner->error) {
	case   0: break;
	case  -1: eslog::error("MKL PDSS: input inconsistent.\n"); break;
	case  -2: eslog::info("MKL PDSS: not enough memory.\n"); ret = false; break;
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

	return ret;
#endif
}

