
#include "transientfirstorderimplicitsolver.h"

#include "esinfo/stepinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"

#include "physics/system/linearsystem.h"
#include "physics/system/builder/builder.h"
#include "physics/substepsolver/substepsolver.h"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

#include "config/ecf/physics/physicssolver/transientfirstorderimplicit.h"
#include "config/ecf/physics/physicssolver/loadstep.h"

#include "math/matrix.h"
#include "math/vector.dense.h"

#include <cmath>

using namespace espreso;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(LinearSystem *system, SubStepSolver *subStepSolver, TransientFirstOrderImplicitSolverConfiguration &configuration, double duration)
: LoadStepSolver(system, subStepSolver, duration),
  _configuration(configuration),
  _alpha(0), _nTimeShift(_configuration.time_step),
  U(NULL), dU(NULL), V(NULL), X(NULL), Y(NULL), dTK(NULL), dTM(NULL)
{
	if (configuration.time_step < 1e-7) {
		eslog::globalerror("Set time step for TRANSIENT solver greater than 1e-7.\n");
	}
}

TransientFirstOrderImplicit::~TransientFirstOrderImplicit()
{
	if (U   != NULL) delete U;
	if (dU  != NULL) delete dU;
	if (V   != NULL) delete V;
	if (X   != NULL) delete X;
	if (Y   != NULL) delete Y;
	if (dTK != NULL) delete dTK;
	if (dTM != NULL) delete dTM;
}

void TransientFirstOrderImplicit::init(LoadStepSolver *previous)
{
	dU = system->solver()->x->shallowCopyStructure();
	U = system->solver()->x->shallowCopyStructure();
	V = system->solver()->x->shallowCopyStructure();
	X = system->solver()->f->shallowCopyStructure();
	Y = system->solver()->f->shallowCopyStructure();
	dTK = system->solver()->x->shallowCopyStructure();
	dTM = system->solver()->x->shallowCopyStructure();

	U->fillData(system->solver()->x);
	if (dynamic_cast<TransientFirstOrderImplicit*>(previous)) {
		V->fillData(dynamic_cast<TransientFirstOrderImplicit*>(previous)->V);
	} else {
		V->fill(0);
	}

	switch (_configuration.method) {
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::CRANK_NICOLSON:
		_alpha = 0.5;
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::GALERKIN:
		_alpha = 2 / 3;
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::BACKWARD_DIFF:
		_alpha = 1;
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::USER:
		_alpha = _configuration.alpha;
		if (_alpha <= 0 || _alpha > 1) {
			eslog::globalerror("Alpha has to be from interval (0, 1>.\n");
		}
		break;
	default:
		eslog::globalerror("Not supported first order implicit solver method.\n");
	}
}

void TransientFirstOrderImplicit::updateStructuralMatrices()
{
	system->builder->matrices &= Builder::Request::K | Builder::Request::M | Builder::Request::RBCf;
	system->assemble();

	if (system->builder->matrices & (Builder::Request::M | Builder::Request::f)) {
		X->sum(1 / (_alpha * step::time.shift), U, (1 - _alpha) / _alpha, V);
		system->assembler()->M->apply(X, Y);
		system->solver()->f->add(1, Y);
		eslog::checkpointln("PHYSICS SOLVER: MATRICES POST-PROCESSED");
	}
}

bool TransientFirstOrderImplicit::runNextSubstep()
{
	step::time.previous = step::time.current;
	step::time.current += _nTimeShift;
	if (step::time.current + step::time.precision >= step::time.final) {
		step::time.current = step::time.final;
	}
	step::time.shift = step::time.current - step::time.previous;
	system->nextSubstep();

	system->builder->internalForceReduction = 1;
	system->builder->timeIntegrationConstantK = 1;
	system->builder->timeIntegrationConstantC = 0;
	system->builder->timeIntegrationConstantM = 1 / (_alpha * step::time.shift);

	switch (_configuration.method) {
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::CRANK_NICOLSON:
		eslog::solver("\n = ==================================== CRANK  NICOLSON ==================================== =\n");
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::GALERKIN:
		eslog::solver("\n = ======================================== GALERKIN ======================================= =\n");
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::BACKWARD_DIFF:
		eslog::solver("\n = ===================================== BACKWARD DIFF ===================================== =\n");
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::USER:
		eslog::solver("\n = =================================== TRANSIENT  SOLVER =================================== =\n");
		break;
	default:
		eslog::globalerror("Not supported first order implicit solver method.\n");
	}
	eslog::solver(" =  LOAD STEP %2d, SUBSTEP %4d, TIME %10.6f, TIME STEP %10.6f, FINAL TIME %10.6f =\n", step::step.loadstep + 1, step::step.substep + 1, step::time.current, step::time.shift, step::time.final);
	eslog::solver(" = ----------------------------------------------------------------------------------------- =\n");

	bool ret = subStepSolver->solve(*this);
	if (!ret) return false;

	dU->sum(1, system->solver()->x, -1, U);
	_nTimeShift = step::time.shift;

	if (_configuration.auto_time_stepping.allowed && step::time.current < step::time.final) {
		eslog::solver(" = -------------------------------- AUTOMATIC TIME STEPPING -------------------------------- =\n");
		if (dU->at(0)->norm() / U->at(0)->norm() < 1e-5) {
			_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * step::time.shift);
		} else {
			system->assembler()->K->apply(dU, dTK);
			system->assembler()->M->apply(dU, dTM);

			double resFreq = dU->at(0)->dot(dTK->at(0)) / dU->at(0)->dot(dTM->at(0));
			double oscilationLimit = step::time.shift * resFreq;
			double t1 = _configuration.auto_time_stepping.oscilation_limit / resFreq;

			if (step::time.shift != t1) {
				if (step::time.shift < t1) {
					if (_configuration.auto_time_stepping.IDFactor * step::time.shift < t1) {
						_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * step::time.shift);
					}
				} else {
					if (step::time.shift / _configuration.auto_time_stepping.IDFactor > t1) {
						_nTimeShift = std::max(_configuration.auto_time_stepping.min_time_step, step::time.shift / _configuration.auto_time_stepping.IDFactor);
					}
				}
			}
			eslog::solver(" = RESPONSE FREQUENCY, OSCILATION LIMIT                             %.5e, %.5e =\n", resFreq, oscilationLimit);
		}

		if (std::fabs(step::time.shift - _nTimeShift) / step::time.shift < step::time.precision) {
			eslog::solver(" = TIME STEP UNCHANGED                                                                       =\n");
		} else {

			if (step::time.shift < _nTimeShift) {
				eslog::solver(" = TIME STEP INCREASED TO NEW VALUE                                                 %8.6f =\n", _nTimeShift);
			} else {
				eslog::solver(" = TIME STEP DECREASED TO NEW VALUE                                                 %8.6f =\n", _nTimeShift);
			}
		}
		eslog::checkpointln("PHYSICS SOLVER: TIME STEP COMPUTED");
	}

	if (step::time.shift - step::time.precision < _nTimeShift) {
		V->sum(1 / (_alpha * step::time.shift), dU, -(1 - _alpha) / _alpha, V);
		U->fillData(system->solver()->x);
		system->processSolution();
	} else {
		system->solver()->x->fillData(U);
		system->solutionChanged();
		step::time.current = step::time.previous;
		--step::step.substep;
	}
	eslog::solver(" = ========================================================================================= =\n");
	eslog::solver(" = ================================================================= run time %12.3f s =\n\n", eslog::duration());

	return true;
}


