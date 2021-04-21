/*
 * LinearSolver.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */
//#include <Driver/DissectionSolver.hpp>
#include "FETISystemSolver.h"
#include "feti/dataholder.h"
#include "physics/system/fetisystem.h"
#include "feti/specific/itersolvers.h"
#include "feti/specific/superclusters.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "autoopt/optimizer.h"
#include "mesh/mesh.h"
#include "mesh/store/domainstore.h"
#include "config/ecf/linearsolver/feti.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "timeeval.h"
#include "wrappers/nvtx/w.nvtx.h"

namespace espreso {

struct FETIDataHolder {
	SuperCluster *cluster;
	IterSolver   *solver;
	DataHolder   holder;

	FETIDataHolder(): cluster(NULL), solver(NULL) {}
	~FETIDataHolder(){
		if (cluster) { delete cluster; }
		if (solver) { delete solver; }
	}
};
}

using namespace espreso;

//FETI::FETI()
//{
//	_inner = new FETIDataHolder();
//	timeEvalMain = new TimeEval("ESPRESO Solver Overall Timing");
//}

FETISystemSolver::FETISystemSolver(FETIConfiguration &configuration, FETISolverData &data)
: configuration(configuration), _data(data), _inner(NULL)
{
	timeEvalMain = new TimeEval("ESPRESO Solver Overall Timing");

	if (configuration.auto_optimization.algorithm != AutoOptimizationConfiguration::ALGORITHM::NONE)
	{
		std::vector<ECFParameter*> opt_parameters;
		opt_parameters = {
			configuration.ecfdescription->getParameter(&configuration.preconditioner),
			configuration.ecfdescription->getParameter(&configuration.iterative_solver),
			configuration.ecfdescription->getParameter(&configuration.regularization),
			configuration.ecfdescription->getParameter(&configuration.redundant_lagrange),
			configuration.ecfdescription->getParameter(&configuration.B0_type),
			configuration.ecfdescription->getParameter(&configuration.scaling),
			// this->configuration.getParameter(&this->configuration.use_schur_complement)
			configuration.ecfdescription->getParameter(&configuration.method)
		};
		this->optimizer = new EvolutionaryOptimizer(configuration.auto_optimization, opt_parameters);
	} else {
		this->optimizer = new EmptyOptimizer();
	}
}

void FETISystemSolver::init()
{
	_inner = new FETIDataHolder();
}

void FETISystemSolver::update()
{
	insertK(configuration, _data.K, _data.origK, _data.N1, _data.N2, _data.RegMat);
	insertB1(_data.B1Dirichlet, _data.B1c, _data.B1Gluing, _data.B1duplication, _data.B1Inequality, _data.B1gap, _data.B1Map);
	insertB0(_data.B0);
	insertRHS(_data.f);

	while(!optimizer->set([&]() {
		if (configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS && 
			configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		{ return false; }
		// Intel MKL ERROR: Parameter 5 was incorrect on entry to MKL_DCSRMV.
		else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
			configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
			configuration.method == FETIConfiguration::METHOD::TOTAL_FETI)
		{ return false; }
		// Intel MKL ERROR: Parameter 5 was incorrect on entry to MKL_DCSRMV.
		else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
			configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
			configuration.method == FETIConfiguration::METHOD::TOTAL_FETI)
		{ return false; }
		// Intel MKL ERROR: Parameter 5 was incorrect on entry to MKL_DCSRMV.
		else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
			configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
			configuration.method == FETIConfiguration::METHOD::TOTAL_FETI)
		{ return false; }
		// Intel MKL ERROR: Parameter 5 was incorrect on entry to MKL_DCSRMV.
		else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
			configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
			configuration.method == FETIConfiguration::METHOD::TOTAL_FETI)
		{ return false; }
		// Intel MKL ERROR: Parameter 5 was incorrect on entry to MKL_DCSRMV.
		else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
			configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
			configuration.method == FETIConfiguration::METHOD::TOTAL_FETI)
		{ return false; }
		// if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::QPCE &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::GMRES &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::GMRES &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::GMRES &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::GMRES &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::GMRES &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::PCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::pipePCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::LUMPED &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::pipePCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::pipePCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::pipePCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		// else if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET &&
		// 	configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::pipePCG &&
		// 	configuration.B0_type == FETIConfiguration::B0_TYPE::KERNELS &&
		// 	configuration.method == FETIConfiguration::METHOD::HYBRID_FETI)
		// { return false; }
		
		int ret = update(configuration);
		if (ret >= 0) return true;

		switch (ret) {
			case -3:
				eslog::info("FETI update: MKL Sparse Solver - Error during solution.\n");
				break;
			default:
				eslog::error("FETI update: Unknown error!\n");
				break;
		}

		return false;
	}));
}

void FETISystemSolver::solve()
{
	solve(configuration, _data.x, _data.y);
	_data.x.averageDuplications();
}

double& FETISystemSolver::precision()
{
	return configuration.precision;
}


FETISystemSolver::~FETISystemSolver() {

	if (_inner && _inner->solver) {
		_inner->solver->preproc_timing.printStatsMPI();
		_inner->solver->timing.printStatsMPI();
		_inner->solver->postproc_timing.printStatsMPI();
		_inner->solver->timeEvalAppa.printStatsMPI();
		_inner->solver->timeEvalProj.printStatsMPI();

		if ( _inner->solver->USE_PREC != FETIConfiguration::PRECONDITIONER::NONE ) {
			_inner->solver->timeEvalPrec.printStatsMPI();
		}
	}

	//TODO: Fix timing:  if ( cluster->USE_HFETI == 1 ) cluster->ShowTiming();

	timeEvalMain->totalTime.endWithBarrier();
	timeEvalMain->printStatsMPI();

	delete timeEvalMain;

	if (_inner) {
		delete _inner;
	}
	if (optimizer) {
		delete optimizer;
	}
}

void FETISystemSolver::insertK(FETIConfiguration &configuration, const MatrixCSRFETI &K, const MatrixCSRFETI &origK, const MatrixDenseFETI &N1, const MatrixDenseFETI &N2, const MatrixCSRFETI &RegMat)
{
	auto setType = [] (SparseMatrix &m, MatrixType type) {
		m.mtype = type;
		switch (type) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			m.type = 'S'; break;
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			m.type = 'S'; break;
		case MatrixType::REAL_UNSYMMETRIC:
			m.type = 'G'; break;
		}
	};

	_inner->holder.decomposition = &K;

	_inner->holder.K.resize(K.domains);
	_inner->holder.origK.resize(K.domains);
	_inner->holder.N1.resize(K.domains);
	_inner->holder.N2.resize(K.domains);
	_inner->holder.origKN1.resize(K.domains);
	_inner->holder.origKN2.resize(K.domains);
	_inner->holder.RegMat.resize(K.domains);
	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		_inner->holder.K[d].rows = K[d].nrows;
		_inner->holder.K[d].cols = K[d].ncols;
		_inner->holder.K[d].nnz = K[d].nnz;
		setType(_inner->holder.K[d], K[d].type);

		_inner->holder.K[d].CSR_I_row_indices.assign(K[d].rows, K[d].rows + K[d].nrows + 1);
		_inner->holder.K[d].CSR_J_col_indices.assign(K[d].cols, K[d].cols + K[d].nnz);
		_inner->holder.K[d].CSR_V_values.assign(K[d].vals, K[d].vals + K[d].nnz);

		if (
				configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
				configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {

			_inner->holder.origK[d].rows = origK[d].nrows;
			_inner->holder.origK[d].cols = origK[d].ncols;
			_inner->holder.origK[d].nnz = origK[d].nnz;
			setType(_inner->holder.origK[d], origK[d].type);

			_inner->holder.origK[d].CSR_I_row_indices.assign(origK[d].rows, origK[d].rows + origK[d].nrows + 1);
			_inner->holder.origK[d].CSR_J_col_indices.assign(origK[d].cols, origK[d].cols + origK[d].nnz);
			_inner->holder.origK[d].CSR_V_values.assign(origK[d].vals, origK[d].vals + origK[d].nnz);
		}

		if (configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
			_inner->holder.N1[d].rows = N1[d].nrows;
			_inner->holder.N1[d].cols = N1[d].ncols;
			_inner->holder.N1[d].nnz = N1[d].nrows * N1[d].ncols;
			_inner->holder.N2[d].rows = N2[d].nrows;
			_inner->holder.N2[d].cols = N2[d].ncols;
			_inner->holder.N2[d].nnz = N2[d].nrows * N2[d].ncols;
			_inner->holder.RegMat[d].rows = RegMat[d].nrows;
			_inner->holder.RegMat[d].cols = RegMat[d].ncols;
			_inner->holder.RegMat[d].nnz = RegMat[d].nnz;

			// from ROWMAYOR to COLMAYOR
			_inner->holder.N1[d].dense_values.clear();
			_inner->holder.N1[d].dense_values.reserve(N1[d].nrows * N1[d].ncols);
			for (esint c = 0; c < N1[d].ncols; ++c) {
				for (esint r = 0; r < N1[d].nrows; ++r) {
					_inner->holder.N1[d].dense_values.push_back(N1[d][r][c]);
				}
			}
			_inner->holder.N2[d].dense_values.clear();
			_inner->holder.N2[d].dense_values.reserve(N2[d].nrows * N2[d].ncols);
			for (esint c = 0; c < N2[d].ncols; ++c) {
				for (esint r = 0; r < N2[d].nrows; ++r) {
					_inner->holder.N2[d].dense_values.push_back(N2[d][r][c]);
				}
			}

			if ((esint)_inner->holder.RegMat[d].CSR_V_values.size() != RegMat[d].nnz) {
				_inner->holder.RegMat[d].CSR_I_row_indices.assign(RegMat[d].rows, RegMat[d].rows + RegMat[d].nrows + 1);
				_inner->holder.RegMat[d].CSR_J_col_indices.assign(RegMat[d].cols, RegMat[d].cols + RegMat[d].nnz);
			}
			_inner->holder.RegMat[d].CSR_V_values.assign(RegMat[d].vals, RegMat[d].vals + RegMat[d].nnz);

			if (RegMat[d].nnz) {
				_inner->holder.K[d].MatAddInPlace(_inner->holder.RegMat[d], 'N', 1);
			}
			_inner->holder.RegMat[d].ConvertToCOO(1);

			// WARN: only non-empty matrix can have type, otherwise solver fails
			if (_inner->holder.N1[d].rows && _inner->holder.N1[d].cols) {
				setType(_inner->holder.N1[d], N1[d].type);
			}
			if (_inner->holder.N2[d].rows && _inner->holder.N2[d].cols) {
				setType(_inner->holder.N2[d], N2[d].type);
			}
			if (_inner->holder.RegMat[d].rows && _inner->holder.RegMat[d].cols) {
				setType(_inner->holder.RegMat[d], RegMat[d].type);
			}

			if (
					configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
					configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {

				_inner->holder.origKN1[d] = _inner->holder.N1[d];
				_inner->holder.origKN2[d] = _inner->holder.N2[d];
			}
		}
	}
}

void FETISystemSolver::insertB1(const MatrixIJVFETI &B1Dirichlet, const VectorDenseFETI &c, const MatrixIJVFETI &B1Gluing, const VectorDenseFETI &duplication, const MatrixIJVFETI &B1Inequality, const VectorDenseFETI &gap, const std::vector<esint> &B1Map)
{
	_inner->holder.B1.resize(B1Dirichlet.domains);
	_inner->holder.B1c.resize(B1Dirichlet.domains);
	_inner->holder.B1duplication.resize(B1Dirichlet.domains);
	_inner->holder.LB.resize(B1Dirichlet.domains);
	#pragma omp parallel for
	for (esint d = 0; d < B1Dirichlet.domains; ++d) {
		_inner->holder.B1[d].rows = B1Inequality[d].nrows;
		_inner->holder.B1[d].cols = B1Dirichlet[d].ncols;
		_inner->holder.B1[d].nnz = B1Dirichlet[d].nnz + B1Gluing[d].nnz + B1Inequality[d].nnz;
		_inner->holder.B1[d].mtype = B1Dirichlet[d].type;

		if ((esint)_inner->holder.B1[d].I_row_indices.size() != B1Dirichlet[d].nnz + B1Gluing[d].nnz) {
			_inner->holder.B1[d].I_row_indices.clear();
			_inner->holder.B1[d].I_row_indices.insert(_inner->holder.B1[d].I_row_indices.end(), B1Dirichlet[d].rows, B1Dirichlet[d].rows + B1Dirichlet[d].nnz);
			_inner->holder.B1[d].I_row_indices.insert(_inner->holder.B1[d].I_row_indices.end(), B1Gluing[d].rows, B1Gluing[d].rows + B1Gluing[d].nnz);
			_inner->holder.B1[d].I_row_indices.insert(_inner->holder.B1[d].I_row_indices.end(), B1Inequality[d].rows, B1Inequality[d].rows + B1Inequality[d].nnz);

			_inner->holder.B1[d].J_col_indices.clear();
			_inner->holder.B1[d].J_col_indices.insert(_inner->holder.B1[d].J_col_indices.end(), B1Dirichlet[d].cols, B1Dirichlet[d].cols + B1Dirichlet[d].nnz);
			_inner->holder.B1[d].J_col_indices.insert(_inner->holder.B1[d].J_col_indices.end(), B1Gluing[d].cols, B1Gluing[d].cols + B1Gluing[d].nnz);
			_inner->holder.B1[d].J_col_indices.insert(_inner->holder.B1[d].J_col_indices.end(), B1Inequality[d].cols, B1Inequality[d].cols + B1Inequality[d].nnz);
		}

		_inner->holder.B1[d].V_values.clear();
		_inner->holder.B1[d].V_values.insert(_inner->holder.B1[d].V_values.end(), B1Dirichlet[d].vals, B1Dirichlet[d].vals + B1Dirichlet[d].nnz);
		_inner->holder.B1[d].V_values.insert(_inner->holder.B1[d].V_values.end(), B1Gluing[d].vals, B1Gluing[d].vals + B1Gluing[d].nnz);
		_inner->holder.B1[d].V_values.insert(_inner->holder.B1[d].V_values.end(), B1Inequality[d].vals, B1Inequality[d].vals + B1Inequality[d].nnz);

		_inner->holder.B1c[d].clear();
		_inner->holder.B1c[d].insert(_inner->holder.B1c[d].end(), c[d].vals, c[d].vals + c[d].size);
		_inner->holder.B1c[d].resize(c[d].size + duplication[d].size, 0);
		_inner->holder.B1c[d].insert(_inner->holder.B1c[d].end(), gap[d].vals, gap[d].vals + gap[d].size);

		_inner->holder.B1duplication[d].clear();
		_inner->holder.B1duplication[d].resize(c[d].size, 1); // dirichlet put always 1
		_inner->holder.B1duplication[d].insert(_inner->holder.B1duplication[d].end(), duplication[d].vals, duplication[d].vals + duplication[d].size);
		_inner->holder.B1duplication[d].resize(_inner->holder.B1duplication[d].size() + gap[d].size, 1); // can it be arbitrary ??

		_inner->holder.LB[d].resize(c[d].size + duplication[d].size, -std::numeric_limits<double>::infinity());
		_inner->holder.LB[d].resize(_inner->holder.LB[d].size() + B1Inequality[d].nnz, 0);
	}

	_inner->holder.B1Map = B1Map;
//	_inner->holder.B1clustersMap.clear();
//	for (esint d = 0; d < B1Dirichlet.domains; ++d) {
//		for (esint i = 0; i < B1Dirichlet[d].nnz; i++) {
//			_inner->holder.B1clustersMap.push_back({ B1Dirichlet[d].rows[i] - 1, info::mpi::rank });
//		}
//		std::sort(_inner->holder.B1clustersMap.begin(), _inner->holder.B1clustersMap.end());
//	}
//
//	for (size_t i = 0; i < B1Map.size(); i += 3) {
//		_inner->holder.B1clustersMap.push_back({ B1Map[i], info::mpi::rank });
//		if (B1Map[i + 1] != info::mpi::rank) {
//			_inner->holder.B1clustersMap.back().push_back(B1Map[i + 1]);
//		}
//		if (B1Map[i + 2] != info::mpi::rank) {
//			_inner->holder.B1clustersMap.back().push_back(B1Map[i + 2]);
//		}
//	}
//	for (esint d = 0; d < B1Dirichlet.domains; ++d) {
//		size_t size = _inner->holder.B1clustersMap.size();
//		for (esint i = 0; i < B1Inequality[d].nnz; i++) {
//			_inner->holder.B1clustersMap.push_back({ B1Inequality[d].rows[i] - 1, info::mpi::rank });
//		}
//		std::sort(_inner->holder.B1clustersMap.begin() + size, _inner->holder.B1clustersMap.end());
//	}
}

void FETISystemSolver::insertB0(const MatrixIJVFETI &B0)
{
	_inner->holder.B0.resize(B0.domains);

	#pragma omp parallel for
	for (esint d = 0; d < B0.domains; ++d) {
		_inner->holder.B0[d].rows = B0[d].nrows;
		_inner->holder.B0[d].cols = B0[d].ncols;
		_inner->holder.B0[d].nnz = B0[d].nnz;
		_inner->holder.B0[d].mtype = B0[d].type;
		switch (_inner->holder.B0[d].mtype) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			_inner->holder.B0[d].type = 'S'; break;
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			_inner->holder.B0[d].type = 'S'; break;
		case MatrixType::REAL_UNSYMMETRIC:
			_inner->holder.B0[d].type = 'G'; break;
		}

		if ((esint)_inner->holder.B0[d].I_row_indices.size() != B0[d].nnz) {
			_inner->holder.B0[d].I_row_indices.assign(B0[d].rows, B0[d].rows + B0[d].nnz);
			_inner->holder.B0[d].J_col_indices.assign(B0[d].cols, B0[d].cols + B0[d].nnz);
		}

		_inner->holder.B0[d].V_values.clear();
		_inner->holder.B0[d].V_values.assign(B0[d].vals, B0[d].vals + B0[d].nnz);
	}
}

void FETISystemSolver::insertRHS(const VectorsDenseFETI &f)
{
	_inner->holder.F.resize(f.holder()->domains);
	#pragma omp parallel for
	for (esint d = 0; d < f.holder()->domains; ++d) {
		_inner->holder.F[d].assign(f[0][d].vals, f[0][d].vals + f[0][d].size);
	}
}

// make partial initialization according to updated matrices
int FETISystemSolver::update(FETIConfiguration &configuration)
{
	// TODO update appropriate solver objects and stop steeling matrices! :)
	// factorization and preconditioners and HFETI preprocessing

	delete _inner->cluster;
	delete _inner->solver;

	//instance->computeKernels(configuration.regularization, configuration.sc_size);

	_inner->holder.B0.resize(info::mesh->domains->size);
	_inner->holder.N1.resize(info::mesh->domains->size);
	_inner->holder.N2.resize(info::mesh->domains->size);
	_inner->holder.RegMat.resize(info::mesh->domains->size);

	std::string type;
	switch (configuration.method) {
	case FETIConfiguration::METHOD::TOTAL_FETI:
		type = "FETI, ";
		break;
	case FETIConfiguration::METHOD::HYBRID_FETI:
		type = "HFETI, ";
		break;
	}
	switch (configuration.iterative_solver) {
	case FETIConfiguration::ITERATIVE_SOLVER::PCG:
		type += "PCG, ";
		break;
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG:
		type += "ORTL PCG, ";
		break;
	case FETIConfiguration::ITERATIVE_SOLVER::pipePCG:
		type += "PIPE PCG, ";
		break;
	case FETIConfiguration::ITERATIVE_SOLVER::GMRES:
		type += "GMRES, ";
		break;
	case FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB:
		type += "BICGSTAB, ";
		break;
	case FETIConfiguration::ITERATIVE_SOLVER::QPCE:
		type += "QPCE, ";
		break;
	case FETIConfiguration::ITERATIVE_SOLVER::PCG_CP:
		type += "PCG CP, ";
		break;
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP:
		type += "ORT PCG CP, ";
		break;
	}
	switch (configuration.preconditioner) {
	case FETIConfiguration::PRECONDITIONER::NONE:
		type += "NONE";
		break;
	case FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION:
		type += "WEIGHTS";
		break;
	case FETIConfiguration::PRECONDITIONER::LUMPED:
		type += "LUMPED";
		break;
	case FETIConfiguration::PRECONDITIONER::DIRICHLET:
		type += "DIRICHLET";
		break;
	case FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET:
		type += "SDIRICHLET";
		break;
	case FETIConfiguration::PRECONDITIONER::MAGIC:
		type += "MAGIC";
		break;
	}

	eslog::solver("     - ---- LINEAR SOLVER -------------------------------------------------------------- -\n");
	eslog::solver("     - | SOLVER :: ESPRESO        TYPE :: %44s | -\n", type.c_str());

	TimeEvent timeSolClusterInit(string("Solver - Cluster init (K factorization incl.)"));timeSolClusterInit.start();
	_inner->cluster = new SuperClusterCPU(configuration, &_inner->holder);
	timeSolClusterInit.endWithBarrier(); timeEvalMain->addEvent(timeSolClusterInit);
	_inner->solver  = new IterSolver(configuration);

	int ret = init(info::mesh->neighbors, configuration);
	if (ret < 0) return ret;

	if (info::ecf->output.print_matrices > 0) {
		eslog::storedata(" STORE MATRICES FOR FETI ITER SOLVER\n");
		std::string prefix = utils::debugDirectory() + "/fetisolver/init";
		_inner->cluster->printInitData(prefix.c_str(), info::ecf->output.print_matrices);
	}

	return ret;
}

// run solver and store primal and dual solution
void FETISystemSolver::solve(FETIConfiguration &configuration, VectorsDenseFETI &x, VectorsDenseFETI &y)
{
	if (
			std::any_of(_inner->holder.K.begin(), _inner->holder.K.end(), [] (const SparseMatrix &K) { return K.mtype == MatrixType::REAL_UNSYMMETRIC; }) &&
			configuration.iterative_solver != FETIConfiguration::ITERATIVE_SOLVER::GMRES &&
			configuration.iterative_solver != FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB) {

		eslog::error("Invalid Linear Solver configuration: Only GMRES and BICGSTAB can solve unsymmetric system.\n");
	}

	double start = eslog::time();
	eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.precision);

	while (!optimizer->run([&] () {
		int ret = Solve(_inner->holder.F, _inner->holder.primalSolution, _inner->holder.dualSolution);
		// Successful run of solver
		if (ret >= 0 && ret < configuration.max_iterations) return true;
		// Solver exceeded the maximum number of iterations and did not converge
		if (ret >= configuration.max_iterations) {
			eslog::info("FETI solve: Maximum number of iterations has been exceeded.\n");
			return false;
		}
		
		// Solver errors
		switch (ret)
		{
			case -1:
				eslog::info("FETI solve: Regular CG with conjugate projector not implemented yet.\n");
				break;
			case -2:
				eslog::info("FETI solve: Geneo requires dirichlet preconditioner.\n");
				break;
			case -3:
				eslog::info("FETI solve: MKL Sparse Solver - Error during solution.\n");
				break;
			default:
				eslog::error("FETI solve: Unknow error!\n");
				break;
		}
		return false;
	}));


	#pragma omp parallel for
	for (size_t d = 0; d < _inner->holder.F.size(); ++d) {
		memcpy(x[0][d].vals, _inner->holder.primalSolution[d].data(), sizeof(double) * _inner->holder.F[d].size());
		memcpy(y[0][d].vals, _inner->holder.dualSolution[d].data(), sizeof(double) * _inner->holder.F[d].size());
	}

	eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
	eslog::solver("     - --------------------------------------------------------------------------------- -\n");
}


int FETISystemSolver::setup_HTFETI() {
// Setup Hybrid FETI part of the solver

	if (_inner->cluster->USE_HFETI == 1) {
		TimeEvent timeHFETIprec(string("Solver - HFETI preprocessing"));
		timeHFETIprec.start();

		int ret = _inner->cluster->SetClusterHFETI();
		if (ret < 0) return ret;

		timeHFETIprec.endWithBarrier();
		timeEvalMain->addEvent(timeHFETIprec);

//		ESLOG(MEMORY) << "After HFETI preprocessing process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	}
	
	return 0;
}

void FETISystemSolver::setup_LocalSchurComplement(FETIConfiguration &configuration) {
// Computation of the Local Schur Complements
//TODO: update for multiple clusters

	if (_inner->cluster->USE_KINV == 1) {
		PUSH_RANGE("Solver - Schur Complement asm.", 1)
		TimeEvent KSCMem(string("Solver - SC asm. mem [MB]")); KSCMem.startWithoutBarrier(GetProcessMemory_u());
		TimeEvent timeSolSC2(string("Solver - Schur Complement asm."));timeSolSC2.start();
		bool USE_FLOAT = false;
		if (configuration.schur_precision == FETIConfiguration::FLOAT_PRECISION::SINGLE) {
			USE_FLOAT = true;
		}
		_inner->cluster->Create_SC_perDomain(USE_FLOAT);
		timeSolSC2.endWithBarrier(); timeEvalMain->addEvent(timeSolSC2);
		KSCMem.endWithoutBarrier(GetProcessMemory_u()); //KSCMem.printLastStatMPIPerNode();
		POP_RANGE // END Solver - Schur Complement asm.
//		 ESLOG(MEMORY) << "After K inv. process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	} else {
		for (size_t d = 0; d < _inner->cluster->domains.size(); d++) {
			_inner->cluster->domains[d]->B1Kplus.is_on_acc = 0;
		}
	}

}

void FETISystemSolver::setup_Preconditioner() {
// Load Matrix K, Regularization
		PUSH_RANGE("Solver - Setup preconditioners", 2)
		TimeEvent timeRegKproc(string("Solver - Setup preconditioners")); timeRegKproc.start();
//		 ESLOG(MEMORY) << "Before - Setup preconditioners - process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent KregMem(string("Solver - Setup preconditioners mem. [MB]")); KregMem.startWithoutBarrier(GetProcessMemory_u());

		_inner->cluster->SetupPreconditioner();

		KregMem.endWithoutBarrier(GetProcessMemory_u()); //KregMem.printLastStatMPIPerNode();

//		 ESLOG(MEMORY) << "After - Setup preconditioners " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
		timeRegKproc.endWithBarrier();
		timeEvalMain->addEvent(timeRegKproc);
		POP_RANGE  // END Solver - Setup preconditioners
}

void FETISystemSolver::setup_FactorizationOfStiffnessMatrices() {
// K Factorization
		PUSH_RANGE("Solver - K factorization", 3)
		TimeEvent timeSolKproc(string("Solver - K factorization")); timeSolKproc.start();
		TimeEvent KFactMem(string("Solver - K factorization mem. [MB]")); KFactMem.startWithoutBarrier(GetProcessMemory_u());

		_inner->cluster->SetupKsolvers();

		KFactMem.endWithoutBarrier(GetProcessMemory_u()); //KFactMem.printLastStatMPIPerNode();
//		 ESLOG(MEMORY) << "After K solver setup process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
		timeSolKproc.endWithBarrier();
		timeEvalMain->addEvent(timeSolKproc);
		POP_RANGE // END Solver - K factorization
}


void FETISystemSolver::setup_SetDirichletBoundaryConditions() {
// Set Dirichlet Boundary Condition
		PUSH_RANGE("Solver - Set Dirichlet Boundary Condition", 4)
		TimeEvent timeSetInitialCondition(string("Solver - Set Dirichlet Boundary Condition"));
		timeSetInitialCondition.start();

		#pragma omp parallel for
		for (size_t d = 0; d < _inner->cluster->domains.size(); d++) {
			_inner->cluster->domains[d]->vec_c  = _inner->holder.B1c[d];
			_inner->cluster->domains[d]->vec_lb = _inner->holder.LB[d];
		}

		timeSetInitialCondition.endWithBarrier();
		timeEvalMain->addEvent(timeSetInitialCondition);
		POP_RANGE // END Solver - Set Dirichlet Boundary Condition
}

void FETISystemSolver::setup_CreateG_GGt_CompressG() {
		PUSH_RANGE("Solver - Create G GGt Compress G", 5)
		TimeEvent timeSolPrec(string("Solver - FETI Preprocessing")); timeSolPrec.start();

//		 ESLOG(MEMORY) << "Solver Preprocessing - HFETI with regularization from K matrix";
//		 ESLOG(MEMORY) << "process " << info::mpi::rank << " uses "<< Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent G1_perCluster_time("Setup G1 per Cluster time   - preprocessing"); G1_perCluster_time.start();

		TimeEvent G1_perCluster_mem("Setup G1 per Cluster memory - preprocessing"); G1_perCluster_mem.startWithoutBarrier(GetProcessMemory_u());
		_inner->cluster->Create_G_perCluster();

		G1_perCluster_time.end(); G1_perCluster_time.printStatMPI();
		G1_perCluster_mem.endWithoutBarrier(GetProcessMemory_u()); G1_perCluster_mem.printStatMPI();

//		 ESLOG(MEMORY) << "Created G1 per cluster";
//		 ESLOG(MEMORY) << "Before HFETI create GGt process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent solver_Preprocessing_time("Setup GGt time   - preprocessing"); solver_Preprocessing_time.start();
		TimeEvent solver_Preprocessing_mem("Setup GGt memory - preprocessing"); solver_Preprocessing_mem.start();

		_inner->solver->Preprocessing(*_inner->cluster);

		solver_Preprocessing_time.end(); solver_Preprocessing_time.printStatMPI();
		solver_Preprocessing_mem.end();  solver_Preprocessing_mem.printStatMPI();

//		 ESLOG(MEMORY) << "Create GGtInv";
//		 ESLOG(MEMORY) << "process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent solver_G1comp_time("Setup G1 compression time   - preprocessing"); solver_G1comp_time.start();
		TimeEvent solver_G1comp_mem("Setup G1 compression memory - preprocessing");  solver_G1comp_mem.start();

		_inner->cluster->Compress_G1(); // Compression of Matrix G1 to work with compressed lambda vectors

		solver_G1comp_time.end(); solver_G1comp_time.printStatMPI();
		solver_G1comp_mem.end();  solver_G1comp_mem.printStatMPI();
//		 ESLOG(MEMORY) << "G1 compression";
//		 ESLOG(MEMORY) << "process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		timeSolPrec.endWithBarrier(); timeEvalMain->addEvent(timeSolPrec);
		POP_RANGE // END Solver - Create G GGt Compress G

}



void FETISystemSolver::setup_InitClusterAndSolver(FETIConfiguration &configuration)
{

	_inner->cluster->SetupCommunicationLayer();


	// *** Iter Solver Set-up ***************************************************************************************************************************
	_inner->solver->CG_max_iter 		= configuration.max_iterations;
	_inner->solver->USE_GGtINV  		= 1;
	_inner->solver->precision 			= configuration.precision;
	_inner->solver->USE_PREC 			= configuration.preconditioner;

	_inner->solver->USE_KINV 			= configuration.use_schur_complement ? 1 : 0;
	_inner->solver->PAR_NUM_THREADS 	= info::env::PAR_NUM_THREADS;
	_inner->solver->SOLVER_NUM_THREADS  = info::env::SOLVER_NUM_THREADS;

	switch (configuration.method) {
	case FETIConfiguration::METHOD::TOTAL_FETI:
		_inner->solver->USE_HFETI = false;
		break;
	case FETIConfiguration::METHOD::HYBRID_FETI:
		_inner->solver->USE_HFETI = true;
		break;
	default:
		eslog::error("Unsupported FETI METHOD.\n");
	}


}

// TODO: const parameters
int FETISystemSolver::init(const std::vector<int> &neighbors, FETIConfiguration &configuration)
{
	//mkl_cbwr_set(MKL_CBWR_COMPATIBLE);

	// Overall Linear Solver Time measurement structure
	timeEvalMain->totalTime.startWithBarrier();



	// *** Initialize the Cluster and Solver structures
	setup_InitClusterAndSolver(configuration);




	// *** Setup B0 matrix
	// setup_B0Matrices();

	// *** Setup B1 matrix
	// setup_B1Matrices();

	// *** Setup R matrix
	// setup_KernelMatrices();

	// *** Set Dirichlet Boundary Condition
	// setup_SetDirichletBoundaryConditions();

	// // *** if TFETI is used or if HTFETI and analytical kernel are used we can compute GGt here - between solution in terms of peak memory
	// if (  !(_inner->cluster->USE_HFETI == 1 && configuration.regularization == FETIConfiguration::REGULARIZATION::ALGEBRAIC)  ) {
	// 	setup_CreateG_GGt_CompressG();
	// }

	// *** Computation of the Schur Complement
	setup_LocalSchurComplement(configuration);

	// *** setup all preconditioners
	setup_Preconditioner();

	// *** K Factorization
	setup_FactorizationOfStiffnessMatrices();

	// *** if KINV is used on GPU, we need LSC first 
	if (  !(_inner->cluster->USE_HFETI == 1 && configuration.regularization == FETIConfiguration::REGULARIZATION::ALGEBRAIC)  ) {
		setup_CreateG_GGt_CompressG();
	}


	// *** Setup Hybrid FETI part of the solver
	int ret = setup_HTFETI();
	if (ret < 0) return ret;

	// *** For algebraic kernels, GGt needs to be computed after HTFETI preprocessing
	if (_inner->cluster->USE_HFETI == 1 && configuration.regularization == FETIConfiguration::REGULARIZATION::ALGEBRAIC) {
		setup_CreateG_GGt_CompressG();
	}

	// Setup Conj Projector

	if ( configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::GENEO ) {
		// solver->CreateConjProjector(*cluster);
	}

	// Cleanup of unnecessary objects
	// TODO: This can be a problem in some cases - need to be verified
	_inner->cluster->my_lamdas_map_indices.clear();

	#pragma omp parallel for
	for (size_t d = 0; d < _inner->cluster->domains.size(); d++) {
		_inner->cluster->domains[d]->B1.Clear();
	}

//	ESLOG(MEMORY) << "End of preprocessing - process " << info::mpi::rank << " uses " << Measure::processMemory() << " MB";
//	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	return 0;
}

int FETISystemSolver::Solve( std::vector < std::vector < double > >  & f_vec,
								std::vector < std::vector < double > >  & prim_solution)
{

	std::vector < std::vector < double > > dual_solution;
	return Solve(f_vec, prim_solution, dual_solution);
}

int FETISystemSolver::Solve( std::vector < std::vector < double > >  & f_vec,
								std::vector < std::vector < double > >  & prim_solution,
								std::vector < std::vector < double > > & dual_solution) {

	PUSH_RANGE("Solver - CG Solver runtime", 6)
	TimeEvent timeSolCG(string("Solver - CG Solver runtime"));
	timeSolCG.start();

	int ret = _inner->solver->Solve    ( *_inner->cluster, f_vec, prim_solution, dual_solution );
	timeSolCG.endWithBarrier();
	timeEvalMain->addEvent(timeSolCG);
	return ret;
	POP_RANGE // END Solver - CG Solver runtime
}

void FETISystemSolver::Postprocessing( ) {

}

void FETISystemSolver::finalize() {

	// Show Linear Solver Runtime Evaluation
//	solver->preproc_timing.printStatsMPI();
//	solver->timing.printStatsMPI();
//	solver->postproc_timing.printStatsMPI();
//	solver->timeEvalAppa.printStatsMPI();
//	solver->timeEvalProj.printStatsMPI();
//
//	if ( solver->USE_PREC != FETIConfiguration::PRECONDITIONER::NONE ) solver->timeEvalPrec.printStatsMPI();
//
//	//TODO: Fix timing:  if ( cluster->USE_HFETI == 1 ) cluster->ShowTiming();
//
//	 timeEvalMain->totalTime.endWithBarrier();
//	 timeEvalMain->printStatsMPI();
}

void FETISystemSolver::CheckSolution( vector < vector < double > > & prim_solution ) {
	// *** Solutin correctnes test **********************************************************************************************

	//	double max_v = 0.0;
//		for (esint i = 0; i < number_of_subdomains_per_cluster; i++)
//			for (size_t j = 0; j < prim_solution[i].size(); j++)
//				if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );
//
//	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.startWithoutBarrier(0.0); max_sol_ev.endWithoutBarrier(max_v);
//
//	double max_vg;
//	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, info::mpi::MPICommunicator );
//	//ESINFO(DETAILS) << "Maxvalue in solution = " << std::setprecision(12) << max_vg;

	//max_sol_ev.printLastStatMPIPerNode(max_vg);
	// *** END - Solutin correctnes test ******************************************************************************************

}

void FETISystemSolver::Preprocessing( std::vector < std::vector < esint > > & lambda_map_sub) {

}
