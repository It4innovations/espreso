
#include "w.pardiso.systemsolver.h"
#include "physics/system/pardisosystem.h"
#include "math/vector.dense.h"
#include "math/vector.dense.distributed.h"
#include "math/matrix.csr.h"
#include "math/matrix.csr.distributed.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"

using namespace espreso;

PARDISOSystemSolver::PARDISOSystemSolver(PARDISOConfiguration &configuration, PARDISOSolverData &data)
: configuration(configuration), _precision(1e-15), _data(data), _K(NULL), _f(NULL), _x(NULL)
{
#ifndef HAVE_PARDISO
	eslog::error("ESPRESO run-time error: cannot call PARDISO library (the library is not linked).\n");
#endif
	if (info::mpi::size > 1) {
		eslog::error("ESPRESO run-time error: only single MPI process is supported by PARDISO library.\n");
	}
}

void PARDISOSystemSolver::init()
{
	_K = new MatrixCSR();
	_f = new VectorDense();
	_x = new VectorDense();
}

void PARDISOSystemSolver::update()
{
	eslog::solver("   - ---- LINEAR SOLVER ------------------------------------------ -\n");
	eslog::solver("   - | SOLVER :: PARDISO                                         | -\n");
	double start = eslog::time();

	bool fillStructure = false;
	if (_K->nrows == 0) {
		fillStructure = true;
	}

	if (_data.K.type != MatrixType::REAL_UNSYMMETRIC) {
		_K->DataMatrixCSR::deepCopy(&_data.K);
		_K->removeLower(_data.K.type);
	} else {
		_K->::Matrix::operator=(_data.K);
		_K->::_DataMatrixCSR::operator=(_data.K);
	}

	switch (_data.K.type) {
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		eslog::solver("   - | MATRIX TYPE ::           REAL SYMMETRIC POSITIVE DEFINITE | -\n"); break;
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
		eslog::solver("   - | MATRIX TYPE ::                  REAL SYMMETRIC INDEFINITE | -\n"); break;
	case MatrixType::REAL_UNSYMMETRIC:
		eslog::solver("   - | MATRIX TYPE ::  REAL NONSYMMETRIC (STRUCTURALY SYMMETRIC) | -\n"); break;
	}

	if (fillStructure) {
		_K->structureUpdated();
		_K->factorizeSymbolic();
		eslog::solver("   - | SYMBOLIC FACTORIZATION                         %8.3f s | -\n", eslog::time() - start);
		start = eslog::time();
	}
	_K->factorizeNumeric();
	eslog::solver("   - | NUMERICAL FACTORIZATION                        %8.3f s | -\n", eslog::time() - start);

	_f->_DataVectorDense::operator=(_data.f[0]);
}

void PARDISOSystemSolver::solve()
{
	double start = eslog::time();
	_x->_DataVectorDense::operator=(_data.x[0]);
	_K->solve(*_f, *_x);
	eslog::solver("   - | SOLVER TIME                                    %8.3f s | -\n", eslog::time() - start);
	eslog::solver("   - ------------------------------------------------------------- -\n");
}

PARDISOSystemSolver::~PARDISOSystemSolver()
{
	if (_K) { delete _K; }
	if (_f) { delete _f; }
	if (_x) { delete _x; }
}

