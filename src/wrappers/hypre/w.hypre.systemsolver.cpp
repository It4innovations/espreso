
#include "w.hypre.systemsolver.h"
#include "physics/system/hypresystem.h"
#include "basis/utilities/communication.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include <vector>
#include <numeric>
#include <fstream>

#ifdef HAVE_HYPRE

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"

namespace espreso {
struct HYPREDataHolder {
	std::vector<esint> rows;
	HYPRE_IJMatrix K;
	HYPRE_IJVector f, x;
	std::vector<HYPRE_IJVector> N;
};
}

#endif

using namespace espreso;

HYPRESystemSolver::HYPRESystemSolver(HYPREConfiguration &configuration, HYPRESolverData &data)
: configuration(configuration), _roffset(0), _nrows(0), _data(data), _inner(NULL)
{
#ifndef HAVE_HYPRE
	eslog::error("ESPRESO run-time error: cannot call HYPRE library (the library is not linked).\n");
#endif
}

void HYPRESystemSolver::init()
{
#ifdef HAVE_HYPRE
	_nrows = _data.K.nrows - _data.K.nhalo;
	_roffset = _nrows;

	_inner = new HYPREDataHolder();
	Communication::exscan(_roffset);

	_inner->rows.resize(_nrows);
	std::iota(_inner->rows.begin(), _inner->rows.end(), _roffset + 1);

	HYPRE_IJMatrixCreate(info::mpi::comm, _roffset + 1, _roffset + _nrows, _roffset + 1, _roffset + _nrows, &_inner->K);
	HYPRE_IJVectorCreate(info::mpi::comm, _roffset + 1, _roffset + _nrows, &_inner->f);
	HYPRE_IJVectorCreate(info::mpi::comm, _roffset + 1, _roffset + _nrows, &_inner->x);

	HYPRE_IJMatrixSetObjectType(_inner->K, HYPRE_PARCSR);
	HYPRE_IJVectorSetObjectType(_inner->f, HYPRE_PARCSR);
	HYPRE_IJVectorSetObjectType(_inner->x, HYPRE_PARCSR);

	HYPRE_IJMatrixInitialize(_inner->K);
	HYPRE_IJVectorInitialize(_inner->f);
	HYPRE_IJVectorInitialize(_inner->x);

	_inner->N.resize(_data.N.nvectors);
	for (esint n = 0; n < _data.N.nvectors; ++n) {
		HYPRE_IJVectorCreate(info::mpi::comm, _roffset + 1, _roffset + _nrows, &_inner->N[n]);
		HYPRE_IJVectorSetObjectType(_inner->N[n], HYPRE_PARCSR);
		HYPRE_IJVectorInitialize(_inner->N[n]);
	}
#endif
}

void HYPRESystemSolver::update()
{
#ifdef HAVE_HYPRE
	std::vector<esint> ncols;
	ncols.reserve(_nrows);
	for (esint r = 0; r < _nrows; r++) {
		ncols.push_back(_data.K.rows[r + _data.K.nhalo + 1] - _data.K.rows[r + _data.K.nhalo]);
	}

	HYPRE_IJMatrixSetValues(
			_inner->K, _nrows, ncols.data(),
			_inner->rows.data(),
			_data.K.cols + _data.K.rows[_data.K.nhalo] - 1,
			_data.K.vals + _data.K.rows[_data.K.nhalo] - 1);
	HYPRE_IJMatrixAssemble(_inner->K);

	HYPRE_IJVectorSetValues(_inner->f, _nrows, _inner->rows.data(), _data.f[0].vals + _data.f[0].nhalo);
	HYPRE_IJVectorAssemble(_inner->f);

	for (esint n = 0; n < _data.N.nvectors; ++n) {
		HYPRE_IJVectorSetValues(_inner->N[n], _nrows, _inner->rows.data(), _data.N[n].vals + _data.N[n].nhalo);
		HYPRE_IJVectorAssemble(_inner->N[n]);
	}
#endif
}

double& HYPRESystemSolver::precision()
{
	switch (configuration.solver_type) {
	case HYPREConfiguration::SOLVER_TYPE::BiCGSTAB:
		return configuration.bicgstab.relative_conv_tol;
	case HYPREConfiguration::SOLVER_TYPE::BoomerAMG:
		return configuration.boomeramg.convergence_tolerance;
	case HYPREConfiguration::SOLVER_TYPE::CGNR:
		return configuration.cgnr.relative_conv_tol;
	case HYPREConfiguration::SOLVER_TYPE::FlexGMRES:
		return configuration.flexgmres.relative_conv_tol;
	case HYPREConfiguration::SOLVER_TYPE::GMRES:
		return configuration.gmres.relative_conv_tol;
	case HYPREConfiguration::SOLVER_TYPE::LGMRES:
		return configuration.lgmres.relative_conv_tol;
	case HYPREConfiguration::SOLVER_TYPE::PCG:
		return configuration.pcg.relative_conv_tol;
	default:
		eslog::globalerror("Required precision of unknown HYPRE solver type.\n");
		exit(0);
	}
}

HYPRESystemSolver::~HYPRESystemSolver()
{
#ifdef HAVE_HYPRE
	HYPRE_IJMatrixDestroy(_inner->K);
	HYPRE_IJVectorDestroy(_inner->f);
	HYPRE_IJVectorDestroy(_inner->x);
	for (size_t n = 0; n < _inner->N.size(); ++n) {
		HYPRE_IJVectorDestroy(_inner->N[n]);
	}
	delete _inner;
#endif
}

#ifdef HAVE_HYPRE
static void setBoomerAMG(HYPRE_Solver &boomerAMG, HYPRESolverData &data, HYPREDataHolder *inner, const HYPREBoomerAMGConfiguration &configuration)
{
	HYPRE_BoomerAMGCreate(&boomerAMG);

	HYPRE_BoomerAMGSetMinIter(boomerAMG,configuration.min_iterations);
	HYPRE_BoomerAMGSetMaxCoarseSize(boomerAMG,configuration.max_coarest_grid_size);
	HYPRE_BoomerAMGSetMinCoarseSize(boomerAMG,configuration.min_coarest_grid_size);
	HYPRE_BoomerAMGSetMaxLevels(boomerAMG,configuration.max_multigrid_levels);

	switch (configuration.coarsening_type) {
		case HYPREBoomerAMGConfiguration::COARSENING_TYPE::CLJP:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 0);
			break;
		case HYPREBoomerAMGConfiguration::COARSENING_TYPE::Ruge_Stuben:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 3);
			break;
		case HYPREBoomerAMGConfiguration::COARSENING_TYPE::Falgout:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 6);
			break;
		case HYPREBoomerAMGConfiguration::COARSENING_TYPE::PMIS:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 8);
			break;
		case HYPREBoomerAMGConfiguration::COARSENING_TYPE::HMIS:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 10);
			break;
		case HYPREBoomerAMGConfiguration::COARSENING_TYPE::CGC:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 21);
			break;
		case HYPREBoomerAMGConfiguration::COARSENING_TYPE::CGC_E:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 22);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}

	switch (configuration.interpolation_type) {
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::CLASSIC_MODIFF:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 0);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::LS:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 1);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::CLASSIC_HYPERBOLIC:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 2);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::DIRECT:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 3);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::MULTLIPASS:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 4);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::MULTIPASSS_SEPARATION:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 5);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::EXTENDED_I:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 6);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::EXTENDED_I_NO_NEIGHBOR:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 7);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::STANDARD:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 8);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::STANDARD_SEPPARATION:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 9);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::CLASSIC_BLOCK:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 10);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::CLASSIC_BLOCK_DIAG:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 11);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::FF:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 12);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::FF1:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 13);
			break;
		case HYPREBoomerAMGConfiguration::INTERPOLATION_TYPE::EXTENDED:
			HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 14);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}

	HYPRE_BoomerAMGSetTruncFactor(boomerAMG, configuration.interp_trunc_factor);
	HYPRE_BoomerAMGSetPMaxElmts(boomerAMG, configuration.max_element_per_row);
	HYPRE_BoomerAMGSetSepWeight(boomerAMG, configuration.weight_separation);

	switch (configuration.cycle_type) {
		case HYPREBoomerAMGConfiguration::CYCLE_TYPE::V_CYCLE:
			HYPRE_BoomerAMGSetCycleType(boomerAMG, 1);
			break;
		case HYPREBoomerAMGConfiguration::CYCLE_TYPE::W_CYCLE:
			HYPRE_BoomerAMGSetCycleType(boomerAMG, 2);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}

	HYPRE_BoomerAMGSetNumSweeps(boomerAMG, configuration.sweeps_num);

	auto getRelaxType = [] (HYPREBoomerAMGConfiguration::RELAX_TYPE type) {
		switch (type) {
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::Jacobi:
				return 0;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::GSS:
				return 1;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::GSIP:
				return 2;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::HGSF:
				return 3;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::HGSB:
				return 4;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::HCHGS:
				return 5;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::HSGS:
				return 6;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::LSHSGS:
				return 8;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::GE:
				return 9;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::CG:
				return 15;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::Chebyshev:
				return 16;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::FCF:
				return 17;
				break;
			case HYPREBoomerAMGConfiguration::RELAX_TYPE::LSJ:
				return 18;
				break;
			default:
				eslog::internalFailure("not implemented interface to the required solver options.\n");
				return -1;
		}
	};

	HYPRE_BoomerAMGSetRelaxType(boomerAMG, getRelaxType(configuration.relax_type));

	if (configuration.allow_cycle_relax_type){
		auto getCycleTypeNum = [] (HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleRelaxTypeConfiguration::RELAX_TYPE_CYCLE type) {
			switch (type) {
				case HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleRelaxTypeConfiguration::RELAX_TYPE_CYCLE::DOWN:
					return 1;
					break;
				case HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleRelaxTypeConfiguration::RELAX_TYPE_CYCLE::UP:
					return 2;
					break;
				case HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleRelaxTypeConfiguration::RELAX_TYPE_CYCLE::COAREST:
					return 3;
					break;
				default:
					eslog::internalFailure("not implemented interface to the required solver options.\n");
					return -1;
			}
		};

		HYPRE_BoomerAMGSetCycleRelaxType(boomerAMG, getRelaxType(configuration.cycle_relax_type.relax_type), getCycleTypeNum(configuration.cycle_relax_type.relax_type_cycle));

	}

	HYPRE_BoomerAMGSetRelaxOrder(boomerAMG, configuration.relax_order);
	HYPRE_BoomerAMGSetRelaxWt(boomerAMG, configuration.relax_weight);

	if (configuration.allow_cycle_num_sweeps){
		auto getCycleTypeNum = [] (HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleSweepsConfiguration::RELAX_TYPE_CYCLE type) {
			switch (type) {
				case HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleSweepsConfiguration::RELAX_TYPE_CYCLE::DOWN:
					return 1;
					break;
				case HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleSweepsConfiguration::RELAX_TYPE_CYCLE::UP:
					return 2;
					break;
				case HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleSweepsConfiguration::RELAX_TYPE_CYCLE::COAREST:
					return 3;
					break;
				default:
					eslog::internalFailure("not implemented interface to the required solver options.\n");
					return -1;
			}
		};

		HYPRE_BoomerAMGSetCycleRelaxType(boomerAMG, configuration.cycle_sweep_spec.sweeps_num_specific, getCycleTypeNum(configuration.cycle_sweep_spec.relax_type_cycle));

	}

	if( getRelaxType(configuration.relax_type)>2 && getRelaxType(configuration.relax_type)<8 ){
	//**************************
	//**************************
		//HYPRE_BoomerAMGSetLevelRelaxWt(boomerAMG, configuration.level_relax_weight,0);
		HYPRE_BoomerAMGSetOuterWt(boomerAMG, configuration.out_relax_weight);
		//HYPRE_BoomerAMGSetLevelOuterWt(boomerAMG, configuration.outer_weight_level,0);
	//**************************
	//**************************
	}

	if ( getRelaxType(configuration.relax_type)==16){
		HYPRE_BoomerAMGSetChebyOrder(boomerAMG,configuration.cheby_order);
		HYPRE_BoomerAMGSetChebyFraction(boomerAMG,configuration.cheby_fraction);
	}

	switch (configuration.smooth_type) {
		case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::SCHWARZ:
			HYPRE_BoomerAMGSetSmoothType(boomerAMG, 6);
			break;
		case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::PILUT:
			HYPRE_BoomerAMGSetSmoothType(boomerAMG, 7);
			break;
		case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::PARASAILS:
			HYPRE_BoomerAMGSetSmoothType(boomerAMG, 8);
			break;
		case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::EUCLID:
			HYPRE_BoomerAMGSetSmoothType(boomerAMG, 9);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}


	auto getSmoothType = [] (HYPREBoomerAMGConfiguration::SMOOTH_TYPE type) {
		switch (type) {
			case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::SCHWARZ:
				return 6;
				break;
			case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::PILUT:
				return 7;
				break;
			case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::PARASAILS:
				return 8;
				break;
			case HYPREBoomerAMGConfiguration::SMOOTH_TYPE::EUCLID:
				return 9;
				break;
			default:
				eslog::internalFailure("not implemented interface to the required solver options.\n");
				return -1;
		}
	};

	HYPRE_BoomerAMGSetSmoothType(boomerAMG, getSmoothType(configuration.smooth_type));

	HYPRE_BoomerAMGSetSmoothNumLevels(boomerAMG, configuration.smooth_level_num);
	HYPRE_BoomerAMGSetSmoothNumSweeps(boomerAMG, configuration.smooth_sweep_num);

	if (getSmoothType(configuration.smooth_type)==6){
		HYPRE_BoomerAMGSetVariant(boomerAMG, configuration.schwarz_variant);

		HYPRE_BoomerAMGSetOverlap(boomerAMG, configuration.schwarz_overlap);
		HYPRE_BoomerAMGSetDomainType(boomerAMG, configuration.schwarz_doamin);
		HYPRE_BoomerAMGSetSchwarzUseNonSymm(boomerAMG, configuration.schwarz_use_nonsym);
	}

	if(configuration.old_default){
		HYPRE_BoomerAMGSetOldDefault(boomerAMG);
	}

	HYPRE_BoomerAMGSetStrongThreshold(boomerAMG, configuration.amg_strength_treshold);
	HYPRE_BoomerAMGSetMaxRowSum(boomerAMG,configuration.diag_dominant_strength);
	HYPRE_BoomerAMGSetSCommPkgSwitch(boomerAMG,configuration.comm_strength_treshold);
	HYPRE_BoomerAMGSetNonGalerkinTol(boomerAMG,configuration.non_galerkin_drop_tol);
	//HYPRE_BoomerAMGSetLevelNonGalerkinTol(boomerAMG,configuration.level_non_galerkin_drop_tol);

	switch (configuration.measure_type) {
		case HYPREBoomerAMGConfiguration::MEASURE_TYPE::LOCAL:
			HYPRE_BoomerAMGSetMeasureType(boomerAMG, 0);
			break;
		case HYPREBoomerAMGConfiguration::MEASURE_TYPE::GLOBAL:
			HYPRE_BoomerAMGSetMeasureType(boomerAMG, 1);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}

	switch (configuration.nodal_sys_coarsening) {
		case HYPREBoomerAMGConfiguration::NODAL_SYS_COARSENING::UNKNOWN:
			HYPRE_BoomerAMGSetNodal(boomerAMG, 0);
			break;
		case HYPREBoomerAMGConfiguration::NODAL_SYS_COARSENING::FROBENIUSNORM:
			HYPRE_BoomerAMGSetNodal(boomerAMG, 1);
			break;
		case HYPREBoomerAMGConfiguration::NODAL_SYS_COARSENING::SUM_ABSOLUTE_BLOCK:
			HYPRE_BoomerAMGSetNodal(boomerAMG, 2);
			break;
		case HYPREBoomerAMGConfiguration::NODAL_SYS_COARSENING::LARGEST_ELEM_BLOCK:
			HYPRE_BoomerAMGSetNodal(boomerAMG, 3);
			break;
		case HYPREBoomerAMGConfiguration::NODAL_SYS_COARSENING::SUM_ROW:
			HYPRE_BoomerAMGSetNodal(boomerAMG, 4);
			break;
		case HYPREBoomerAMGConfiguration::NODAL_SYS_COARSENING::ALL_BLOC_SUM:
			HYPRE_BoomerAMGSetNodal(boomerAMG, 6);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}

	HYPRE_BoomerAMGSetNodalDiag(boomerAMG, configuration.nodal_diag_treatment);
	HYPRE_BoomerAMGSetCGCIts(boomerAMG, configuration.cgc_pathes);

	HYPRE_BoomerAMGSetAggNumLevels(boomerAMG, configuration.aggressive_coarsening_levels);
	HYPRE_BoomerAMGSetNumPaths(boomerAMG, configuration.aggressive_coarsening_degree);
	HYPRE_BoomerAMGSetAggTruncFactor(boomerAMG, configuration.aggressive_coarsening_trunc_factor);
	HYPRE_BoomerAMGSetAggP12TruncFactor(boomerAMG, configuration.aggressive_coarsening_p12_trunc_fact);
	HYPRE_BoomerAMGSetAggPMaxElmts(boomerAMG, configuration.aggressive_coarsening_p_max_elements);
	HYPRE_BoomerAMGSetAggP12MaxElmts(boomerAMG, configuration.aggressive_coarsening_p12_max_elements);

	switch (configuration.aggressive_coarsening_interp) {
		case HYPREBoomerAMGConfiguration::AGGRESSIVE_COARSENING_INTERP::TWO_STAGE_IEXTENDED:
			HYPRE_BoomerAMGSetAggInterpType(boomerAMG, 1);
			break;
		case HYPREBoomerAMGConfiguration::AGGRESSIVE_COARSENING_INTERP::TWO_STAGE_STANDARD:
			HYPRE_BoomerAMGSetAggInterpType(boomerAMG, 2);
			break;
		case HYPREBoomerAMGConfiguration::AGGRESSIVE_COARSENING_INTERP::TWO_STAGE_EXTENDED:
			HYPRE_BoomerAMGSetAggInterpType(boomerAMG, 3);
			break;
		case HYPREBoomerAMGConfiguration::AGGRESSIVE_COARSENING_INTERP::MULTIPASS:
			HYPRE_BoomerAMGSetAggInterpType(boomerAMG, 4);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}

	HYPRE_BoomerAMGSetNumFunctions(boomerAMG, data.numFnc);
	switch (configuration.interp_vec_variant) {
	case HYPREBoomerAMGConfiguration::INTERP_VEC_VARIANT::NONE:
		HYPRE_BoomerAMGSetInterpVecVariant(boomerAMG, 0); break;
	case HYPREBoomerAMGConfiguration::INTERP_VEC_VARIANT::GM1:
		HYPRE_BoomerAMGSetInterpVecVariant(boomerAMG, 1); break;
	case HYPREBoomerAMGConfiguration::INTERP_VEC_VARIANT::GM2:
		HYPRE_BoomerAMGSetInterpVecVariant(boomerAMG, 2); break;
	case HYPREBoomerAMGConfiguration::INTERP_VEC_VARIANT::LM:
		HYPRE_BoomerAMGSetInterpVecVariant(boomerAMG, 3); break;
	}
	HYPRE_BoomerAMGSetInterpVecQMax(boomerAMG, 0);
	HYPRE_BoomerAMGSetInterpVecAbsQTrunc(boomerAMG, 0);
}

static void setBoomerAMGSolver(HYPRE_Solver &boomerAMG, HYPRESolverData &data, HYPREDataHolder *inner, const HYPREBoomerAMGConfiguration &configuration)
{
	setBoomerAMG(boomerAMG, data, inner, configuration);
	HYPRE_BoomerAMGSetTol(boomerAMG,configuration.convergence_tolerance);
	HYPRE_BoomerAMGSetMaxIter(boomerAMG,configuration.max_iterations);
}

static void setBoomerAMGPreconditioner(HYPRE_Solver &boomerAMG, HYPRESolverData &data, HYPREDataHolder *inner, const HYPREBoomerAMGConfiguration &configuration)
{
	setBoomerAMG(boomerAMG, data, inner, configuration);
	HYPRE_BoomerAMGSetTol(boomerAMG, 0);
	HYPRE_BoomerAMGSetMaxIter(boomerAMG,1);
}

static void setParaSailsPreconditioner(HYPRE_Solver &parasails, const HYPREParaSailsConfiguration &configuration)
{
	HYPRE_ParaSailsSetParams(parasails, configuration.threshold, configuration.n_levels);
	HYPRE_ParaSailsSetFilter(parasails, configuration.filter);

	switch (configuration.symmetry) {
		case HYPREParaSailsConfiguration::SYMMETRY::NON_INF:
			 HYPRE_ParaSailsSetSym(parasails, 0);
			break;
		case HYPREParaSailsConfiguration::SYMMETRY::SPD:
			HYPRE_ParaSailsSetSym(parasails, 1);
			break;
		case HYPREParaSailsConfiguration::SYMMETRY::NON_DEF_SPD:
			HYPRE_ParaSailsSetSym(parasails, 2);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
	}

	HYPRE_ParaSailsSetLoadbal(parasails, configuration.loadbal);
	HYPRE_ParaSailsSetReuse(parasails, configuration.reuse);
	HYPRE_ParaSailsSetLogging(parasails, configuration.logging);
}


static void setEuclidPreconditioner(HYPRE_Solver &euclid, const HYPREEuclidConfiguration &configuration)
{
	HYPRE_EuclidSetLevel(euclid,configuration.levels);
	HYPRE_EuclidSetBJ(euclid,configuration.set_bj);
	HYPRE_EuclidSetStats(euclid,configuration.stats);
	HYPRE_EuclidSetMem(euclid,configuration.memory_stats);
	HYPRE_EuclidSetSparseA(euclid,configuration.sparse_tol);
	HYPRE_EuclidSetRowScale(euclid,configuration.row_scale);
	HYPRE_EuclidSetILUT(euclid,configuration.ilut_tol);
}

static void setPilutPreconditioner(HYPRE_Solver &pilut, const HYPREPilutConfiguration &configuration)
{
	HYPRE_ParCSRPilutSetMaxIter(pilut, configuration.max_iter);
	HYPRE_ParCSRPilutSetDropTolerance(pilut, configuration.drop_tol);
	HYPRE_ParCSRPilutSetFactorRowSize(pilut, configuration.row_size);
}
#endif

void HYPRESystemSolver::solve()
{
#ifdef HAVE_HYPRE
	std::fill(_data.x[0].vals + _data.x[0].nhalo, _data.x[0].vals + _data.x[0].size, 0);
	HYPRE_IJVectorSetValues(_inner->x, _nrows, _inner->rows.data(), _data.x[0].vals + _data.x[0].nhalo);
	HYPRE_IJVectorAssemble(_inner->x);

	eslog::solver("     - ---- LINEAR SOLVER -------------------------------------------------------------- -\n");
	double start = eslog::time();

	esint iterations;
	double norm;

	HYPRE_ParCSRMatrix K;
	HYPRE_ParVector f, x;
	std::vector<HYPRE_ParVector> N(_data.N.nvectors);
	HYPRE_IJMatrixGetObject(_inner->K, (void**) &K);
	HYPRE_IJVectorGetObject(_inner->f, (void**) &f);
	HYPRE_IJVectorGetObject(_inner->x, (void**) &x);

	HYPRE_Solver solver;
	HYPRE_Solver preconditioner;
	switch (configuration.solver_type) {
	case HYPREConfiguration::SOLVER_TYPE::BoomerAMG:

		setBoomerAMGSolver(solver, _data, _inner, configuration.boomeramg);
		HYPRE_BoomerAMGSetLogging(solver, 0);
		switch (configuration.boomeramg.solver_info) {
		case HYPREBoomerAMGConfiguration::SOLVER_INFO::NO_INFO:
			HYPRE_BoomerAMGSetPrintLevel(solver, 0);
			break;
		case HYPREBoomerAMGConfiguration::SOLVER_INFO::SETUP_INFO:
			HYPRE_BoomerAMGSetPrintLevel(solver, 1);
			break;
		case HYPREBoomerAMGConfiguration::SOLVER_INFO::SOLVE_INFO:
			HYPRE_BoomerAMGSetPrintLevel(solver, 2);
			break;
		case HYPREBoomerAMGConfiguration::SOLVER_INFO::SETUP_SOLVE_INFO:
			HYPRE_BoomerAMGSetPrintLevel(solver, 3);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
		}

		eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                                BOOMERAMG | -\n");
		eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.boomeramg.convergence_tolerance);
		HYPRE_BoomerAMGSetup(solver, K, f, x);
		HYPRE_BoomerAMGSolve(solver, K, f, x);

		HYPRE_BoomerAMGGetNumIterations(solver, &iterations);
		HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &norm);

		eslog::solver("     - | SOLUTION RESIDUAL NORM                                           %e | -\n", norm);
		eslog::solver("     - | ITERATIONS TO CONVERGENCE                                            %8d | -\n", iterations);
		if (norm > configuration.boomeramg.convergence_tolerance) {
			eslog::warning("     - |                       >>> SOLVER DOES NOT CONVERGED <<<                       | -\n");
		}
		eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
		eslog::solver("     - --------------------------------------------------------------------------------- -\n");

		HYPRE_BoomerAMGDestroy(solver);

		break;
	case HYPREConfiguration::SOLVER_TYPE::PCG:
		HYPRE_ParCSRPCGCreate(info::mpi::comm, &solver);

		HYPRE_PCGSetMaxIter(solver, configuration.pcg.max_iterations);
		HYPRE_PCGSetTol(solver, configuration.pcg.relative_conv_tol);
		HYPRE_PCGSetAbsoluteTol(solver, configuration.pcg.absolute_conv_tol);
		HYPRE_PCGSetResidualTol(solver, configuration.pcg.residual_conv_tol);
		HYPRE_PCGSetTwoNorm(solver, configuration.pcg.two_norm?1:0);
		HYPRE_PCGSetRecomputeResidual(solver, configuration.pcg.recompute_residual_end?1:0);
		HYPRE_PCGSetRecomputeResidualP(solver, configuration.pcg.recompute_residual_p?1:0);

		HYPRE_PCGSetLogging(solver, 0);

		switch (configuration.pcg.solver_info) {
		case HYPREPCGConfiguration::SOLVER_INFO::NO_INFO:
			HYPRE_PCGSetPrintLevel(solver, 0);
			break;
		case HYPREPCGConfiguration::SOLVER_INFO::SETUP_INFO:
			HYPRE_PCGSetPrintLevel(solver, 1);
			break;
		case HYPREPCGConfiguration::SOLVER_INFO::SOLVE_INFO:
			HYPRE_PCGSetPrintLevel(solver, 2);
			break;
		case HYPREPCGConfiguration::SOLVER_INFO::SETUP_SOLVE_INFO:
			HYPRE_PCGSetPrintLevel(solver, 3);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
		}

		switch (configuration.pcg.preconditioner) {
		case HYPREPCGConfiguration::PRECONDITIONER::BoomerAMG:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                           PCG, BOOMERAMG | -\n");
			setBoomerAMGPreconditioner(preconditioner, _data, _inner, configuration.pcg.boomeramg);
			if (_data.N.nvectors) {
				for (esint n = 0; n < _data.N.nvectors; ++n) {
					HYPRE_IJVectorGetObject(_inner->N[n], (void**)(N.data() + n));
				}
				HYPRE_BoomerAMGSetInterpVectors(preconditioner, _data.N.nvectors, N.data());
			}
			HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::ParaSails:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                           PCG, PARASAILS | -\n");
			HYPRE_ParaSailsCreate(info::mpi::comm, &preconditioner);
			setParaSailsPreconditioner(preconditioner, configuration.pcg.parasails);
			HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::Euclid:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                              PCG, EUCLID | -\n");
			HYPRE_EuclidCreate(info::mpi::comm, &preconditioner);
			setEuclidPreconditioner(preconditioner, configuration.pcg.euclid);
			HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::Pilut:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                               PCG, PILUT | -\n");
			HYPRE_ParCSRPilutCreate(info::mpi::comm, &preconditioner);
			setPilutPreconditioner(preconditioner, configuration.pcg.pilut);
			HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup, preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}

		eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.pcg.relative_conv_tol);
		HYPRE_ParCSRPCGSetup(solver, K, f, x);
		HYPRE_ParCSRPCGSolve(solver, K, f, x);

		HYPRE_PCGGetNumIterations(solver, &iterations);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &norm);
		eslog::solver("     - | SOLUTION RESIDUAL NORM                                           %e | -\n", norm);
		eslog::solver("     - | ITERATIONS TO CONVERGENCE                                            %8d | -\n", iterations);
		if (norm > configuration.pcg.relative_conv_tol) {
			eslog::warning("     - |                       >>> SOLVER DOES NOT CONVERGED <<<                       | -\n");
		}
		eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
		eslog::solver("     - --------------------------------------------------------------------------------- -\n");

		HYPRE_ParCSRPCGDestroy(solver);

		switch (configuration.pcg.preconditioner) {
		case HYPREPCGConfiguration::PRECONDITIONER::BoomerAMG:
			HYPRE_BoomerAMGDestroy(preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::ParaSails:
			HYPRE_ParaSailsDestroy(preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::Euclid:
			HYPRE_EuclidDestroy(preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::Pilut:
			HYPRE_ParCSRPilutDestroy(preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}
		break;

	case HYPREConfiguration::SOLVER_TYPE::GMRES:
		HYPRE_ParCSRGMRESCreate(info::mpi::comm, &solver);

		HYPRE_GMRESSetMaxIter(solver, configuration.gmres.max_iterations);
		HYPRE_GMRESSetTol(solver, configuration.gmres.relative_conv_tol);
		HYPRE_GMRESSetAbsoluteTol(solver, configuration.gmres.absolute_conv_tol);
		HYPRE_GMRESSetKDim(solver, configuration.gmres.restarts);

		HYPRE_GMRESSetLogging(solver, 0);

		switch (configuration.gmres.solver_info) {
		case HYPREGMRESConfiguration::SOLVER_INFO::NO_INFO:
			HYPRE_GMRESSetPrintLevel(solver, 0);
			break;
		case HYPREGMRESConfiguration::SOLVER_INFO::SETUP_INFO:
			HYPRE_GMRESSetPrintLevel(solver, 1);
			break;
		case HYPREGMRESConfiguration::SOLVER_INFO::SOLVE_INFO:
			HYPRE_GMRESSetPrintLevel(solver, 2);
			break;
		case HYPREGMRESConfiguration::SOLVER_INFO::SETUP_SOLVE_INFO:
			HYPRE_GMRESSetPrintLevel(solver, 3);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
		}

		switch (configuration.gmres.preconditioner) {
		case HYPREGMRESConfiguration::PRECONDITIONER::BoomerAMG:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                         GMRES, BOOMERAMG | -\n");
			setBoomerAMGPreconditioner(preconditioner, _data, _inner, configuration.gmres.boomeramg);
			HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::ParaSails:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                         GMRES, PARASAILS | -\n");
			HYPRE_ParaSailsCreate(info::mpi::comm, &preconditioner);
			setParaSailsPreconditioner(preconditioner, configuration.gmres.parasails);
			HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::Euclid:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                            GMRES, EUCLID | -\n");
			HYPRE_EuclidCreate(info::mpi::comm, &preconditioner);
			setEuclidPreconditioner(preconditioner, configuration.gmres.euclid);
			HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::Pilut:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                             GMRES, PILUT | -\n");
			HYPRE_ParCSRPilutCreate(info::mpi::comm, &preconditioner);
			setPilutPreconditioner(preconditioner, configuration.pcg.pilut);
			HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup, preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::NONE:
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}

		eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.gmres.relative_conv_tol);
		HYPRE_ParCSRGMRESSetup(solver, K, f, x);
		HYPRE_ParCSRGMRESSolve(solver, K, f, x);

		HYPRE_GMRESGetNumIterations(solver, &iterations);
		HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &norm);
		eslog::solver("     - | SOLUTION RESIDUAL NORM                                           %e | -\n", norm);
		eslog::solver("     - | ITERATIONS TO CONVERGENCE                                            %8d | -\n", iterations);
		if (norm > configuration.gmres.relative_conv_tol) {
			eslog::warning("     - |                       >>> SOLVER DOES NOT CONVERGED <<<                       | -\n");
		}
		eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
		eslog::solver("     - --------------------------------------------------------------------------------- -\n");

		HYPRE_ParCSRGMRESDestroy(solver);
		switch (configuration.gmres.preconditioner) {
		case HYPREGMRESConfiguration::PRECONDITIONER::BoomerAMG:
			HYPRE_BoomerAMGDestroy(preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::ParaSails:
			HYPRE_ParaSailsDestroy(preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::Euclid:
			HYPRE_EuclidDestroy(preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::Pilut:
			HYPRE_ParCSRPilutDestroy(preconditioner);
			break;
		case HYPREGMRESConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}
		break;

	case HYPREConfiguration::SOLVER_TYPE::FlexGMRES:
		HYPRE_ParCSRFlexGMRESCreate(info::mpi::comm, &solver);

		HYPRE_FlexGMRESSetMaxIter(solver, configuration.flexgmres.max_iterations);
		HYPRE_FlexGMRESSetTol(solver, configuration.flexgmres.relative_conv_tol);
		HYPRE_FlexGMRESSetAbsoluteTol(solver, configuration.flexgmres.absolute_conv_tol);
		HYPRE_FlexGMRESSetKDim(solver, configuration.flexgmres.restarts);

		HYPRE_FlexGMRESSetLogging(solver, 0);

		switch (configuration.flexgmres.solver_info) {
		case HYPREFlexGMRESConfiguration::SOLVER_INFO::NO_INFO:
			HYPRE_FlexGMRESSetPrintLevel(solver, 0);
			break;
		case HYPREFlexGMRESConfiguration::SOLVER_INFO::SETUP_INFO:
			HYPRE_FlexGMRESSetPrintLevel(solver, 1);
			break;
		case HYPREFlexGMRESConfiguration::SOLVER_INFO::SOLVE_INFO:
			HYPRE_FlexGMRESSetPrintLevel(solver, 2);
			break;
		case HYPREFlexGMRESConfiguration::SOLVER_INFO::SETUP_SOLVE_INFO:
			HYPRE_FlexGMRESSetPrintLevel(solver, 3);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
		}

		switch (configuration.flexgmres.preconditioner) {
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::BoomerAMG:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                     FLEXGMRES, BOOMERAMG | -\n");
			setBoomerAMGPreconditioner(preconditioner, _data, _inner, configuration.flexgmres.boomeramg);
			HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
			break;
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::ParaSails:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                     FLEXGMRES, PARASAILS | -\n");
			HYPRE_ParaSailsCreate(info::mpi::comm, &preconditioner);
			setParaSailsPreconditioner(preconditioner, configuration.flexgmres.parasails);
			HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, preconditioner);
			break;
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::Euclid:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                        FLEXGMRES, EUCLID | -\n");
			HYPRE_EuclidCreate(info::mpi::comm, &preconditioner);
			setEuclidPreconditioner(preconditioner, configuration.flexgmres.euclid);
			HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, preconditioner);
			break;		
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::Pilut:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                         FLEXGMRES, PILUT | -\n");
			HYPRE_ParCSRPilutCreate(info::mpi::comm, &preconditioner);
			setPilutPreconditioner(preconditioner, configuration.pcg.pilut);
			HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup, preconditioner);
			break;
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}

		eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.flexgmres.relative_conv_tol);
		HYPRE_ParCSRFlexGMRESSetup(solver, K, f, x);
		HYPRE_ParCSRFlexGMRESSolve(solver, K, f, x);

		HYPRE_FlexGMRESGetNumIterations(solver, &iterations);
		HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &norm);
		eslog::solver("     - | SOLUTION RESIDUAL NORM                                           %e | -\n", norm);
		eslog::solver("     - | ITERATIONS TO CONVERGENCE                                            %8d | -\n", iterations);
		if (norm > configuration.flexgmres.relative_conv_tol) {
			eslog::warning("     - |                       >>> SOLVER DOES NOT CONVERGED <<<                       | -\n");
		}
		eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
		eslog::solver("     - --------------------------------------------------------------------------------- -\n");

		HYPRE_ParCSRFlexGMRESDestroy(solver);
		switch (configuration.flexgmres.preconditioner) {
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::BoomerAMG:
			HYPRE_BoomerAMGDestroy(preconditioner);
			break;
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::ParaSails:
			HYPRE_ParaSailsDestroy(preconditioner);
			break;
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::Euclid:
			HYPRE_EuclidDestroy(preconditioner);
			break;
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::Pilut:
			HYPRE_ParCSRPilutDestroy(preconditioner);
			break;
		case HYPREFlexGMRESConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}
		break;

	case HYPREConfiguration::SOLVER_TYPE::LGMRES:
		HYPRE_ParCSRLGMRESCreate(info::mpi::comm, &solver);

		HYPRE_LGMRESSetMaxIter(solver, configuration.lgmres.max_iterations);
		HYPRE_LGMRESSetTol(solver, configuration.lgmres.relative_conv_tol);
		HYPRE_LGMRESSetAbsoluteTol(solver, configuration.lgmres.absolute_conv_tol);
		HYPRE_LGMRESSetKDim(solver, configuration.lgmres.restarts);

		HYPRE_LGMRESSetLogging(solver, 0);

		switch (configuration.lgmres.solver_info) {
		case HYPRELGMRESConfiguration::SOLVER_INFO::NO_INFO:
			HYPRE_LGMRESSetPrintLevel(solver, 0);
			break;
		case HYPRELGMRESConfiguration::SOLVER_INFO::SETUP_INFO:
			HYPRE_LGMRESSetPrintLevel(solver, 1);
			break;
		case HYPRELGMRESConfiguration::SOLVER_INFO::SOLVE_INFO:
			HYPRE_LGMRESSetPrintLevel(solver, 2);
			break;
		case HYPRELGMRESConfiguration::SOLVER_INFO::SETUP_SOLVE_INFO:
			HYPRE_LGMRESSetPrintLevel(solver, 3);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver options.\n");
		}

		switch (configuration.lgmres.preconditioner) {
		case HYPRELGMRESConfiguration::PRECONDITIONER::BoomerAMG:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                        LGMRES, BOOMERAMG | -\n");
			setBoomerAMGPreconditioner(preconditioner, _data, _inner, configuration.lgmres.boomeramg);
			HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::ParaSails:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                        LGMRES, PARASAILS | -\n");
			HYPRE_ParaSailsCreate(info::mpi::comm, &preconditioner);
			setParaSailsPreconditioner(preconditioner, configuration.lgmres.parasails);
			HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::Euclid:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                           LGMRES, EUCLID | -\n");
			HYPRE_EuclidCreate(info::mpi::comm, &preconditioner);
			setEuclidPreconditioner(preconditioner, configuration.lgmres.euclid);
			HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::Pilut:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                            LGMRES, PILUT | -\n");
			HYPRE_ParCSRPilutCreate(info::mpi::comm, &preconditioner);
			setPilutPreconditioner(preconditioner, configuration.pcg.pilut);
			HYPRE_LGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup, preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::NONE:
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}

		eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.lgmres.relative_conv_tol);
		HYPRE_ParCSRLGMRESSetup(solver, K, f, x);
		HYPRE_ParCSRLGMRESSolve(solver, K, f, x);

		HYPRE_LGMRESGetNumIterations(solver, &iterations);
		HYPRE_LGMRESGetFinalRelativeResidualNorm(solver, &norm);
		eslog::solver("     - | SOLUTION RESIDUAL NORM                                           %e | -\n", norm);
		eslog::solver("     - | ITERATIONS TO CONVERGENCE                                            %8d | -\n", iterations);
		if (norm > configuration.lgmres.relative_conv_tol) {
			eslog::warning("     - |                       >>> SOLVER DOES NOT CONVERGED <<<                       | -\n");
		}
		eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
		eslog::solver("     - --------------------------------------------------------------------------------- -\n");

		HYPRE_ParCSRLGMRESDestroy(solver);
		switch (configuration.lgmres.preconditioner) {
		case HYPRELGMRESConfiguration::PRECONDITIONER::BoomerAMG:
			HYPRE_BoomerAMGDestroy(preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::ParaSails:
			HYPRE_ParaSailsDestroy(preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::Euclid:
			HYPRE_EuclidDestroy(preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::Pilut:
			HYPRE_ParCSRPilutDestroy(preconditioner);
			break;
		case HYPRELGMRESConfiguration::PRECONDITIONER::NONE:
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}
		break;

	case HYPREConfiguration::SOLVER_TYPE::BiCGSTAB:
		HYPRE_ParCSRBiCGSTABCreate(info::mpi::comm, &solver);

		HYPRE_BiCGSTABSetMaxIter(solver, configuration.bicgstab.max_iterations);
		HYPRE_BiCGSTABSetTol(solver, configuration.bicgstab.relative_conv_tol);
		HYPRE_BiCGSTABSetAbsoluteTol(solver, configuration.bicgstab.absolute_conv_tol);

		HYPRE_BiCGSTABSetLogging(solver, 0);

		switch (configuration.bicgstab.solver_info) {
		case HYPREBiCGSTABConfiguration::SOLVER_INFO::NO_INFO:
			HYPRE_BiCGSTABSetPrintLevel(solver, 0);
			break;
		case HYPREBiCGSTABConfiguration::SOLVER_INFO::SETUP_INFO:
			HYPRE_BiCGSTABSetPrintLevel(solver, 1);
			break;
		case HYPREBiCGSTABConfiguration::SOLVER_INFO::SOLVE_INFO:
			HYPRE_BiCGSTABSetPrintLevel(solver, 2);
			break;
		case HYPREBiCGSTABConfiguration::SOLVER_INFO::SETUP_SOLVE_INFO:
			HYPRE_BiCGSTABSetPrintLevel(solver, 3);
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}

		switch (configuration.bicgstab.preconditioner) {
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::BoomerAMG:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                      BICGSTAB, BOOMERAMG | -\n");
			setBoomerAMGPreconditioner(preconditioner, _data, _inner, configuration.bicgstab.boomeramg);
			HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::ParaSails:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                      BICGSTAB, PARASAILS | -\n");
			HYPRE_ParaSailsCreate(info::mpi::comm, &preconditioner);
			setParaSailsPreconditioner(preconditioner, configuration.bicgstab.parasails);
			HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::Euclid:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                         BICGSTAB, EUCLID | -\n");
			HYPRE_EuclidCreate(info::mpi::comm, &preconditioner);
			setEuclidPreconditioner(preconditioner, configuration.bicgstab.euclid);
			HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::Pilut:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                          BICGSTAB, PILUT | -\n");
			HYPRE_ParCSRPilutCreate(info::mpi::comm, &preconditioner);
			setPilutPreconditioner(preconditioner, configuration.pcg.pilut);
			HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup, preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}

		eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.bicgstab.relative_conv_tol);
		HYPRE_ParCSRBiCGSTABSetup(solver, K, f, x);
		HYPRE_ParCSRBiCGSTABSolve(solver, K, f, x);

		HYPRE_BiCGSTABGetNumIterations(solver, &iterations);
		HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &norm);
		eslog::solver("     - | SOLUTION RESIDUAL NORM                                           %e | -\n", norm);
		eslog::solver("     - | ITERATIONS TO CONVERGENCE                                            %8d | -\n", iterations);
		if (norm > configuration.bicgstab.relative_conv_tol) {
			eslog::warning("     - |                       >>> SOLVER DOES NOT CONVERGED <<<                       | -\n");
		}
		eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
		eslog::solver("     - --------------------------------------------------------------------------------- -\n");

		HYPRE_ParCSRBiCGSTABDestroy(solver);
		switch (configuration.bicgstab.preconditioner) {
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::BoomerAMG:
			HYPRE_BoomerAMGDestroy(preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::ParaSails:
			HYPRE_ParaSailsDestroy(preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::Euclid:
			HYPRE_EuclidDestroy(preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::Pilut:
			HYPRE_ParCSRPilutDestroy(preconditioner);
			break;
		case HYPREBiCGSTABConfiguration::PRECONDITIONER::NONE:
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}
		break;

	case HYPREConfiguration::SOLVER_TYPE::CGNR:
		HYPRE_ParCSRCGNRCreate(info::mpi::comm, &solver);

		HYPRE_CGNRSetMaxIter(solver, configuration.cgnr.max_iterations);
		HYPRE_CGNRSetTol(solver, configuration.cgnr.relative_conv_tol);
		HYPRE_CGNRSetLogging(solver, 0);

		switch (configuration.cgnr.preconditioner) {
		case HYPRECGNRConfiguration::PRECONDITIONER::BoomerAMG:
			eslog::solver("     - | SOLVER ::   HYPRE            TYPE ::                          CGNR, BOOMERAMG | -\n");
			setBoomerAMGPreconditioner(preconditioner, _data, _inner, configuration.cgnr.boomeramg);
		//	HYPRE_CGNRSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
			break;
		case HYPRECGNRConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}

		eslog::solver("     - | REQUESTED STOPPING CRITERIA                                      %e | -\n", configuration.cgnr.relative_conv_tol);
		HYPRE_ParCSRCGNRSetup(solver, K, f, x);
		HYPRE_ParCSRCGNRSolve(solver, K, f, x);

		HYPRE_CGNRGetNumIterations(solver, &iterations);
		HYPRE_CGNRGetFinalRelativeResidualNorm(solver, &norm);
		eslog::solver("     - | SOLUTION RESIDUAL NORM                                           %e | -\n", norm);
		eslog::solver("     - | ITERATIONS TO CONVERGENCE                                            %8d | -\n", iterations);
		if (norm > configuration.cgnr.relative_conv_tol) {
			eslog::warning("     - |                       >>> SOLVER DOES NOT CONVERGED <<<                       | -\n");
		}
		eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
		eslog::solver("     - --------------------------------------------------------------------------------- -\n");

		HYPRE_ParCSRCGNRDestroy(solver);
		switch (configuration.cgnr.preconditioner) {
		case HYPRECGNRConfiguration::PRECONDITIONER::BoomerAMG:
			HYPRE_BoomerAMGDestroy(preconditioner);
			break;
		case HYPRECGNRConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::internalFailure("not implemented interface to the required solver.\n");
		}
		break;
	default:
		eslog::internalFailure("not implemented interface to the required solver.\n");
	}

	HYPRE_IJVectorGetValues(_inner->x, _nrows, _inner->rows.data(), _data.x[0].vals + _data.x[0].nhalo);
	_data.x[0].scatterToUpper();
#endif
}
