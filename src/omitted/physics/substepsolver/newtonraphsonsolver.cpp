
#include "newtonraphsonsolver.h"
#include "esinfo/stepinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"

#include "physics/system/linearsystem.h"
#include "physics/system/builder/builder.h"
#include "physics/loadstepsolver/loadstepsolver.h"

#include "basis/containers/serializededata.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/statisticsstore.h"
#include "config/ecf/physics/physicssolver/loadstep.h"
#include "config/ecf/physics/physicssolver/nonlinear.h"

#include "math/vector.sparse.h"
#include "math/matrix.h"

using namespace espreso;

NewtonRaphson::NewtonRaphson(LinearSystem *system, NonLinearSolverConfiguration &configuration)
: SubStepSolver(system), _configuration(configuration),
  K(NULL), U(NULL), R(NULL), f(NULL), BC(NULL), lsSolution(NULL), lsRHS(NULL), lsResidual(NULL)
{

}

void NewtonRaphson::init(SubStepSolver *previous)
{
	if (U == NULL) {
		U = _system->solver()->x->shallowCopyStructure();
	}

	if (_configuration.check_second_residual) {
		if (R == NULL) {
			R = _system->solver()->f->shallowCopyStructure();
		}
		if (f == NULL) {
			f = _system->solver()->f->shallowCopyStructure();
		}
		if (!_system->solver()->hasReactionForces()) {
			if (K == NULL) {
				K = _system->solver()->K->shallowCopyStructure();
			}
		}
		BC = _system->solver()->BC->shallowCopyStructure();
	} else {
		if (K != NULL) { delete K; K = NULL; }
		if (R != NULL) { delete R; R = NULL; }
		if (f != NULL) { delete f; f = NULL; }
	}

	if (_configuration.line_search) {
		if (lsSolution == NULL) {
			lsSolution = _system->solver()->x->shallowCopyStructure();
		}
		if (lsRHS == NULL) {
			lsRHS = _system->solver()->f->shallowCopyStructure();
		}
		if (lsResidual == NULL) {
			lsResidual = _system->solver()->f->shallowCopyStructure();
		}
	}

	if (previous) {
		if (K) { K->fillData(dynamic_cast<NewtonRaphson*>(previous)->K); }
		if (U) { U->fillData(dynamic_cast<NewtonRaphson*>(previous)->U); }
		if (R) { R->fillData(dynamic_cast<NewtonRaphson*>(previous)->R); }
		if (f) { f->fillData(dynamic_cast<NewtonRaphson*>(previous)->f); }
		if (lsSolution) { lsSolution->fillData(dynamic_cast<NewtonRaphson*>(previous)->lsSolution); }
		if (lsRHS) { lsRHS->fillData(dynamic_cast<NewtonRaphson*>(previous)->lsRHS); }
		if (lsResidual) { lsResidual->fillData(dynamic_cast<NewtonRaphson*>(previous)->lsResidual); }
	}
}

NewtonRaphson::~NewtonRaphson()
{
	if (K != NULL) delete K;
	if (U != NULL) delete U;
	if (R != NULL) delete R;
	if (f != NULL) delete f;
	if (BC != NULL) delete BC;
	if (lsSolution != NULL) { delete lsSolution; }
	if (lsRHS != NULL) { delete lsRHS; }
	if (lsResidual != NULL) { delete lsResidual; }
}

bool NewtonRaphson::hasSameMode(const LoadStepSolverConfiguration &configuration) const
{
	return configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR;
}

void NewtonRaphson::solve(LoadStepSolver &loadStepSolver)
{
	if (!_configuration.check_first_residual && !_configuration.check_second_residual) {
		eslog::globalerror("Turn on at least one convergence parameter for NONLINEAR solver.\n");
	}

	switch (_configuration.method) {
	case NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON:
		eslog::solver("     - -------------------------------- NEWTON  RAPHSON -------------------------------- -\n");
		break;
	case NonLinearSolverConfiguration::METHOD::MODIFIED_NEWTON_RAPHSON:
		eslog::solver("     - ---------------------------- MODIFIED NEWTON RAPHSON ---------------------------- -\n");
		break;
	}

	double solutionNorm = 10 * _configuration.requested_first_residual;

	auto computeDual = [&] () {
		if (!_system->solver()->hasReactionForces()) {
			K->apply(_system->solver()->x, _system->solver()->y);
			_system->solver()->y->sum(1, f, -1, _system->solver()->y);
			eslog::checkpointln("PHYSICS SOLVER: DUAL SOLUTION COMPUTED");
		}
	};

	auto lineSearch = [&] () {
		double alpha = 1;

		double a = 0, b = 1;
		double fa = 0, fb = 0, fx = 0, faStart = 0;

		_system->builder->matrices = Builder::Request::R;

		for (size_t i = 0; i < 6; i++) {
			lsSolution->sum(1, U, alpha, _system->solver()->x);

			lsSolution->swap(_system->solver()->x);
			_system->solutionChanged();
			loadStepSolver.updateStructuralMatrices();
			lsSolution->swap(_system->solver()->x);

			lsResidual->sum(1, lsRHS, -1, _system->solver()->R);

			if (i == 0) {
				faStart = _system->solver()->x->at(0)->dot(_system->solver()->f->at(0));
				fb = _system->solver()->x->at(0)->dot(lsResidual->at(0));
				if ((faStart < 0 && fb < 0) || (faStart >= 0 && fb >= 0)) {
					return alpha;
				}
				fa = faStart;
			} else {
				fx = _system->solver()->x->at(0)->dot(lsResidual->at(0));
				if (fa * fx < 0) {
					b = alpha;
					fb = fx;
				} else if (fb * fx < 0) {
					a = alpha;
					fa = fx;
				}

				if (fabs(fx) <= 0.5 * faStart) {
					alpha = a - fa * ((b - a ) / (fb - fa));
					break;
				}
			}

			alpha = a - fa * ((b - a ) / (fb - fa));
		}

		if (alpha < 0.1) {
			alpha = 0.1;
		}
		if (alpha > .99) {
			alpha = 1;
		}

		_system->solver()->x->scale(alpha);
		eslog::checkpointln("PHYSICS SOLVER: LINE SEARCH PROCESSED");
		return alpha;
	};

	step::step.iteration = 0;
	if (step::isInitial()) {
		_system->builder->tangentMatrixCorrection = false;
		_system->builder->matrices = Builder::Request::KCM | Builder::Request::f | Builder::Request::BC;
		loadStepSolver.updateStructuralMatrices();
		if (_configuration.check_second_residual && !_system->solver()->hasReactionForces()) {
			K->fillData(_system->solver()->K);
			f->fillData(_system->solver()->f);
		}
		_system->setDirichlet();

		eslog::solver("     - INITIAL STEP                              REASSEMBLED MATRICES  ::  %c, %c, %c, %c, %c -    \n",
					(_system->builder->matrices & Builder::Request::K)  ? 'K' : ' ',
					(_system->builder->matrices & Builder::Request::M)  ? 'M' : ' ',
					(_system->builder->matrices & Builder::Request::C)  ? 'C' : ' ',
					(_system->builder->matrices & Builder::Request::R)  ? 'R' : ' ',
					(_system->builder->matrices & Builder::Request::f)  ? 'f' : ' ');

		_system->solve();
		if (_configuration.check_second_residual) {
			computeDual();
		}
		_system->solutionChanged();
	}

	_system->builder->tangentMatrixCorrection = _configuration.tangent_matrix_correction;
	while (step::step.iteration++ < _configuration.max_iterations) {
		U->fillData(_system->solver()->x);
		if (_configuration.method == NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON) {
			_system->builder->matrices = Builder::Request::KCM | Builder::Request::RBCf;
		} else {
			_system->builder->matrices = Builder::Request::RBCf;
		}
		loadStepSolver.updateStructuralMatrices();
		if (_configuration.line_search) {
			lsRHS->fillData(_system->solver()->f);
		}
		if (_configuration.check_second_residual && !_system->solver()->hasReactionForces()) {
			K->fillData(_system->solver()->K);
		}
		if (_configuration.check_second_residual) {
			f->fillData(_system->solver()->f);
		}

		_system->solver()->f->add(-1, _system->solver()->R);

		if (_configuration.check_second_residual) {
			R->sum(1, _system->solver()->f, -1, _system->solver()->y);
			double residualNumerator = R->at(0)->norm();
			BC->fillData(_system->solver()->R);
			R->fillData(f);
			R->fillData(BC);
			double residualDenominator = std::max(R->at(0)->norm(), 1e-3);
			if (residualNumerator / residualDenominator < _configuration.requested_second_residual && step::step.iteration > 1 ) {
				eslog::solver("     - HEAT NORM, CRITERIA                           %.5e / %.5e CONVERGED -\n", residualNumerator, residualDenominator * _configuration.requested_second_residual);
				if (_configuration.check_first_residual) {
					if (solutionNorm < _configuration.requested_first_residual) {
						break;
					}
				} else {
					break;
				}
			} else {
				eslog::solver("     - HEAT NORM, CRITERIA                           %.5e / %.5e           -\n", residualNumerator, residualDenominator * _configuration.requested_second_residual);
			}
			f->fillData(_system->solver()->f);
		}

		_system->solver()->BC->add(-1, U);
		_system->setDirichlet();

		if (_configuration.adaptive_precision) {
//			double solutionPrecisionError = 1;
//			if (time::iteration > 1) {
//				solutionPrecisionError = solutionNumerator / solutionDenominator;
//				solutionPrecision = std::min(_configuration.r_tol * solutionPrecisionError, _configuration.c_fact * solutionPrecision);
//			}
//			eslog::solver("   - ADAPTIVE PRECISION, EPS ERROR       %.5e / %.5e -\n", solutionPrecision, solutionPrecisionError);
		}

		eslog::solver("\n     - EQUIL. ITERATION :: %3d                     REASSEMBLED MATRICES :: %c, %c, %c, %c, %c -\n",
				step::step.iteration,
				(_system->builder->matrices & Builder::Request::K) ? 'K' : ' ',
				(_system->builder->matrices & Builder::Request::M) ? 'M' : ' ',
				(_system->builder->matrices & Builder::Request::C) ? 'C' : ' ',
				(_system->builder->matrices & Builder::Request::R) ? 'R' : ' ',
				(_system->builder->matrices & Builder::Request::f) ? 'f' : ' ');
		_system->solve();
		if (_configuration.check_second_residual) {
			computeDual();
		}

		if (_configuration.line_search) {
			double lsInc = _system->solver()->x->at(0)->absmax();
			double lsAlpha = lineSearch();
			eslog::solver("     - LINE SEARCH, MAX DOF INCREMENT                               %.5f, %.5e -\n", lsAlpha, lsInc);
		}

		if (!_configuration.check_first_residual) {
			_system->solver()->x->add(1, U);
			_system->solutionChanged();
		} else {
			double solutionNumerator = _system->solver()->x->at(0)->norm();
			_system->solver()->x->add(1, U);
			_system->solutionChanged();
			double solutionDenominator = std::max(_system->solver()->x->at(0)->norm(), 1e-3);
			solutionNorm = solutionNumerator / solutionDenominator;

			if (solutionNorm > _configuration.requested_first_residual) {
				eslog::solver("     - TEMP NORM, CRITERIA                           %.5e / %.5e           -\n", solutionNumerator, solutionDenominator * _configuration.requested_first_residual);
			} else {
				eslog::solver("     - TEMP NORM, CRITERIA                           %.5e / %.5e CONVERGED -\n", solutionNumerator, solutionDenominator * _configuration.requested_first_residual);
				if (!_configuration.check_second_residual) {
					break;
				}
			}
		}
	}

	eslog::solver("     - --------------------------------------------------------------------------------- -\n");
}


