
#include "w.superlu.systemsolver.h"
#include "wrappers/mpi/communication.h"
#include "config/ecf/linearsolver/superlu.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "physics/system/superlusystem.h"

#include <vector>

#ifdef HAVE_SUPERLU
#include "superlu_ddefs.h"

namespace espreso {
struct SuperLUDataHolder {
	superlu_dist_options_t options;
	SuperLUStat_t stat;
	SuperMatrix A;
	dScalePermstruct_t ScalePermstruct;
	dLUstruct_t LUstruct;
	dSOLVEstruct_t SOLVEstruct;
	gridinfo_t grid;
	double   *berr;
	int m, n;
	int nprow, npcol;
	int iam, info, ldb, ldx, nrhs;
	double *solution;
	double *rhsValues;
	std::vector<esint> rowData, colData;
};
}
#endif

using namespace espreso;

SuperLUSystemSolver::SuperLUSystemSolver(SuperLUConfiguration &configuration, SuperLUSolverData &data)
: configuration(configuration), _roffset(0), _nrows(0), _precision(1e-15), _data(data), _inner(NULL)
{
#ifndef HAVE_SUPERLU
	eslog::globalerror("ESPRESO run-time error: cannot call SuperLU solver (the library with the solver is not linked).\n");
#endif
}

void SuperLUSystemSolver::init()
{
#ifdef HAVE_SUPERLU
	_nrows = _data.K.nrows - _data.K.nhalo;
	_roffset = _nrows;

	_inner = new SuperLUDataHolder();
	_inner->n = Communication::exscan(_roffset);

	_inner->nprow = configuration.np_row;
	_inner->npcol = configuration.np_col;
	if (_inner->nprow * _inner->npcol != info::mpi::size && info::mpi::rank == 0) {
		eslog::globalerror("ESPRESO: SuperLU error - number of processes in SuperLU grid must equal to number of MPI ranks.\n");
	}
	_inner->nrhs = 1;

	_inner->m = _inner->n;
	_inner->ldb = _nrows;

	superlu_gridinit(info::mpi::comm, _inner->nprow, _inner->npcol, &_inner->grid);
	_inner->iam = _inner->grid.iam;

	set_default_options_dist(&_inner->options);

	_inner->options.SymPattern = yes_no_t::YES;

	dScalePermstructInit(_inner->m, _inner->n, &_inner->ScalePermstruct);
	dLUstructInit(_inner->n, &_inner->LUstruct);
#endif
}

bool SuperLUSystemSolver::update()
{
#ifdef HAVE_SUPERLU
	eslog::solver("     - ---- LINEAR SOLVER -------------------------------------------------------------- -\n");
	eslog::solver("     - | SOLVER :: SUPER LU                      TYPE :: PARALLEL DIRECT SPARSE SOLVER | -\n");
	// TODO: print relevant data
	// eslog::solver("     - | MATRIX TYPE ::                               REAL SYMMETRIC POSITIVE DEFINITE | -\n");
	double start = eslog::time();

	const bool fillStructure = !_inner->rowData.size();
	int nnz_loc = _data.K.rows[_data.K.nrows] - _data.K.rows[_data.K.nhalo];

	_inner->colData.clear();
	_inner->colData.reserve(nnz_loc);
	_inner->rowData.clear();
	_inner->rowData.reserve(_nrows + 1);
	_inner->rowData.push_back(0);

	for (esint i = _data.K.nhalo + 1; i <= _data.K.nrows; i++) {
		_inner->rowData.push_back(_data.K.rows[i] - _data.K.rows[_data.K.nhalo]);
	}

	for (esint i = 0; i < nnz_loc; i++) {
		_inner->colData.push_back(_data.K.cols[_data.K.rows[_data.K.nhalo] - 1 + i] - 1);
	}

	if (!fillStructure) {
		_inner->options.Fact = fact_t::SamePattern_SameRowPerm;
	}

	dCreate_CompRowLoc_Matrix_dist(
			&_inner->A, _inner->m, _inner->n,
			nnz_loc, _nrows, _roffset,
			_data.K.vals + _data.K.rows[_data.K.nhalo] - 1, _inner->colData.data(), _inner->rowData.data(),
			SLU_NR_loc, SLU_D, SLU_GE);

	_inner->rhsValues = _data.f[0].vals + _data.f[0].nhalo;

	if (!(_inner->berr = doubleMalloc_dist(_inner->nrhs))) {
		ABORT("SuperLU: Malloc fails for berr[].");
	}

	eslog::solver("     - | PREPROCESSING                                                      %8.3f s | -\n", eslog::time() - start);
#endif
	return true;
}

bool SuperLUSystemSolver::solve()
{
#ifdef HAVE_SUPERLU
	_inner->solution = _data.x[0].vals + _data.x[0].nhalo;
	double start = eslog::time();
	memcpy(_inner->solution, _inner->rhsValues, _nrows * sizeof(double));

	PStatInit(&_inner->stat);


	pdgssvx(&_inner->options, &_inner->A, &_inner->ScalePermstruct,
		_inner->solution, _inner->ldb, _inner->nrhs, &_inner->grid,
		&_inner->LUstruct, &_inner->SOLVEstruct, _inner->berr, &_inner->stat,
		&_inner->info);

	_data.x.scatterToUpper();

	PStatPrint(&_inner->options, &_inner->stat, &_inner->grid);

	PStatFree(&_inner->stat);
	// Destroy_CompRowLoc_Matrix_dist(&_inner->A);
	Destroy_SuperMatrix_Store_dist(&_inner->A);
	dZeroLblocks(_inner->iam, _inner->n, &_inner->grid, &_inner->LUstruct);
	// Destroy_LU(_inner->n, &_inner->grid, &_inner->LUstruct);

	eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
	eslog::solver("     - --------------------------------------------------------------------------------- -\n");

#endif
	return true;
}

SuperLUSystemSolver::~SuperLUSystemSolver()
{
#ifdef HAVE_SUPERLU
    dScalePermstructFree(&_inner->ScalePermstruct);
    dLUstructFree(&_inner->LUstruct);
    if ( _inner->options.SolveInitialized ) {
        dSolveFinalize(&_inner->options, &_inner->SOLVEstruct);
    }
    SUPERLU_FREE(_inner->berr);

    superlu_gridexit(&_inner->grid);
	delete _inner;
#endif
}
