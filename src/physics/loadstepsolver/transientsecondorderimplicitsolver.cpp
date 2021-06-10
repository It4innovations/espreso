
#include "transientsecondorderimplicitsolver.h"

#include "esinfo/stepinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/evaluator/evaluator.h"

#include "physics/system/linearsystem.h"
#include "physics/system/builder/builder.h"
#include "physics/substepsolver/substepsolver.h"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

#include "config/ecf/physics/physicssolver/transientsecondorderimplicit.h"
#include "config/ecf/physics/physicssolver/loadstep.h"

#include "math/matrix.h"
#include "math/vector.dense.h"

#include <cmath>

using namespace espreso;

TransientSecondOrderImplicit::TransientSecondOrderImplicit(LinearSystem *system, SubStepSolver *subStepSolver, TransientSecondOrderImplicitSolverConfiguration &configuration, double duration)
: LoadStepSolver(system, subStepSolver, duration),
  _configuration(configuration),
  _alpha(_configuration.alpha), _delta(_configuration.delta),
  _massDamping(0), _stiffnessDamping(0),
  _nTimeShift(_configuration.time_step),
  U(NULL), dU(NULL), V(NULL), W(NULL), X(NULL), Y(NULL), Z(NULL), dTK(NULL), dTM(NULL)
{
	if (configuration.time_step < 1e-7) {
		eslog::globalerror("Set time step for TRANSIENT solver greater than 1e-7.\n");
	}

	_alpha += _configuration.numerical_damping;
	_delta *= (1 + _configuration.numerical_damping) * (1 + _configuration.numerical_damping);
	updateConstants();
	updateDamping();
}

TransientSecondOrderImplicit::~TransientSecondOrderImplicit()
{
	if (U   != NULL) delete U;
	if (dU  != NULL) delete dU;
	if (V   != NULL) delete V;
	if (W   != NULL) delete W;
	if (X   != NULL) delete X;
	if (Y   != NULL) delete Y;
	if (Z   != NULL) delete Z;
	if (dTK != NULL) delete dTK;
	if (dTM != NULL) delete dTM;
}

void TransientSecondOrderImplicit::init(LoadStepSolver *previous)
{
	U = system->solver()->x->shallowCopyStructure();
	dU = system->solver()->x->shallowCopyStructure();
	V = system->solver()->f->shallowCopyStructure();
	W = system->solver()->f->shallowCopyStructure();
	X = system->solver()->f->shallowCopyStructure();
	Y = system->solver()->f->shallowCopyStructure();
	Z = system->solver()->f->shallowCopyStructure();
	dTK = system->solver()->x->shallowCopyStructure();
	dTM = system->solver()->x->shallowCopyStructure();

	U->fillData(system->solver()->x);
	if (dynamic_cast<TransientSecondOrderImplicit*>(previous)) {
		V->fillData(dynamic_cast<TransientSecondOrderImplicit*>(previous)->V);
		W->fillData(dynamic_cast<TransientSecondOrderImplicit*>(previous)->W);
	} else {
		V->fill(0);
		W->fill(0);
	}
}

void TransientSecondOrderImplicit::updateConstants()
{
	_newmarkConsts[0] = 1. / (_delta * _nTimeShift * _nTimeShift);
	_newmarkConsts[1] = _alpha / (_delta * _nTimeShift);
	_newmarkConsts[2] = 1. / (_delta * _nTimeShift);
	_newmarkConsts[3] = 1. / (2 * _delta) - 1;
	_newmarkConsts[4] = _alpha / _delta - 1;
	_newmarkConsts[5] = _nTimeShift / 2 * (_alpha / _delta - 2);
	_newmarkConsts[6] = _nTimeShift * (1 - _alpha);
	_newmarkConsts[7] = _nTimeShift * _alpha;
}

void TransientSecondOrderImplicit::updateDamping()
{
	switch (_configuration.damping.rayleigh.type) {
	case RayleighDampingConfiguration::TYPE::NONE:
		break;
	case RayleighDampingConfiguration::TYPE::DIRECT:
		_stiffnessDamping = _configuration.damping.rayleigh.direct_damping.stiffness.evaluator->eval({});
		_massDamping = _configuration.damping.rayleigh.direct_damping.mass.evaluator->eval({});
		break;
	case RayleighDampingConfiguration::TYPE::DAMPING_RATIO: {
		double ratio = _configuration.damping.rayleigh.ratio_damping.ratio.evaluator->eval({});
		double frequency = _configuration.damping.rayleigh.ratio_damping.frequency.evaluator->eval({});
		_stiffnessDamping = 2 * ratio * 2 * M_PI * frequency;
		_massDamping = 2 * ratio / (2 * M_PI * frequency);
	} break;
	}
}

void TransientSecondOrderImplicit::updateStructuralMatrices()
{
	system->builder->matrices &= Builder::Request::K | Builder::Request::M | Builder::Request::RBCf;
	system->assemble();

	switch (_configuration.damping.rayleigh.type) {
	case RayleighDampingConfiguration::TYPE::NONE:
		if (system->builder->matrices & (Builder::Request::M | Builder::Request::f)) {
			X->sum(_newmarkConsts[0], U, _newmarkConsts[2], V);
			X->add(_newmarkConsts[3], W);
			system->assembler()->M->apply(X, Y);
			system->solver()->f->add(1, Y);
			eslog::checkpointln("PHYSICS SOLVER: MATRICES POST-PROCESSED");
		}
		break;
	case RayleighDampingConfiguration::TYPE::DIRECT:
	case RayleighDampingConfiguration::TYPE::DAMPING_RATIO:
		if (system->builder->matrices & (Builder::Request::K | Builder::Request::M | Builder::Request::f)) {
			X->sum(_newmarkConsts[0], U, _newmarkConsts[2], V);
			X->add(_newmarkConsts[3], W);
			X->add(_massDamping * _newmarkConsts[1], U);
			X->add(_massDamping * _newmarkConsts[4], V);
			X->add(_massDamping * _newmarkConsts[5], W);
			system->assembler()->M->apply(X, Y);
			system->solver()->f->add(1, Y);

			X->sum(_stiffnessDamping * _newmarkConsts[1], U, _stiffnessDamping * _newmarkConsts[4], V);
			X->add(_stiffnessDamping * _newmarkConsts[5], W);
			system->assembler()->K->apply(X, Y);
			system->solver()->f->add(1, Y);
			eslog::checkpointln("PHYSICS SOLVER: MATRICES POST-PROCESSED");
		}
		break;
	}
}

void TransientSecondOrderImplicit::runNextSubstep()
{
	step::time.previous = step::time.current;
	step::time.current += _nTimeShift;
	if (step::time.current + step::time.precision >= step::time.final) {
		step::time.current = step::time.final;
	}
	step::time.shift = step::time.current - step::time.previous;
	system->nextSubstep();

	switch (_configuration.damping.rayleigh.type) {
	case RayleighDampingConfiguration::TYPE::NONE:
		break;
	case RayleighDampingConfiguration::TYPE::DIRECT:
		if (	_configuration.damping.rayleigh.direct_damping.stiffness.evaluator->isTimeDependent() ||
				_configuration.damping.rayleigh.direct_damping.mass.evaluator->isTimeDependent()) {
			updateDamping();
		}
		break;
	case RayleighDampingConfiguration::TYPE::DAMPING_RATIO:
		if (	_configuration.damping.rayleigh.ratio_damping.ratio.evaluator->isTimeDependent() ||
				_configuration.damping.rayleigh.ratio_damping.frequency.evaluator->isTimeDependent()) {
			updateDamping();
		}
		break;
	}

	system->builder->internalForceReduction = 1;
	system->builder->timeIntegrationConstantK = 1 + _newmarkConsts[1] * _stiffnessDamping;
	system->builder->timeIntegrationConstantC = 0;
	system->builder->timeIntegrationConstantM = _newmarkConsts[0] + _newmarkConsts[1] * _massDamping;

	switch (_configuration.method) {
	case TransientSecondOrderImplicitSolverConfiguration::METHOD::NEWMARK:
		eslog::solver("\n = ======================================== NEWMARK ======================================== =\n");
		break;
	default:
		eslog::globalerror("Not supported first order implicit solver method.\n");
	}
	eslog::solver(" =  LOAD STEP %2d, SUBSTEP %4d, TIME %10.6f, TIME STEP %10.6f, FINAL TIME %10.5f =\n", step::step.loadstep + 1, step::step.substep + 1, step::time.current, step::time.shift, step::time.final);
	eslog::solver(" = ----------------------------------------------------------------------------------------- =\n");

	subStepSolver->solve(*this);

	dU->sum(1, system->solver()->x, -1, U);
	_nTimeShift = step::time.shift;

	bool changeConstants = false;

	if (_configuration.auto_time_stepping.allowed && step::time.current < step::time.final) {
		eslog::solver(" = -------------------------------- AUTOMATIC TIME STEPPING -------------------------------- =\n");
		if (false && dU->at(0)->norm() / U->at(0)->norm() < 1e-5) {
			_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * step::time.shift);
		} else {
			system->assembler()->K->apply(dU, dTK);
			system->assembler()->M->apply(dU, dTM);

			double resFreq = std::sqrt(std::fabs(dU->at(0)->dot(dTK->at(0)) / (4 * M_PI * M_PI * dU->at(0)->dot(dTM->at(0)))));
			double resPeriod = 1 / resFreq;
			double t2 = resPeriod / _configuration.auto_time_stepping.points_per_period;

			if (step::time.shift != t2) {
				if (step::time.shift < t2) {
					if (_configuration.auto_time_stepping.IDFactor * step::time.shift < t2) {
						_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * step::time.shift);
					}
				} else {
					if (step::time.shift / _configuration.auto_time_stepping.IDFactor > t2) {
						_nTimeShift = std::max(_configuration.auto_time_stepping.min_time_step, step::time.shift / _configuration.auto_time_stepping.IDFactor);
					}
				}
			}

			eslog::solver(" = RESPONSE FREQUENCY, POINTS PER PERIOD                                 %.5e, %6.2f =\n", resFreq, resPeriod / step::time.shift);
		}

		if (std::fabs(step::time.shift - _nTimeShift) / step::time.shift < step::time.precision) {
			eslog::solver(" = TIME STEP UNCHANGED                                                                       =\n");
		} else {
			if (step::time.shift - step::time.precision < _nTimeShift) {
				eslog::solver(" = TIME STEP INCREASED TO NEW VALUE                                                 %8.6f =\n", _nTimeShift);
			} else {
				eslog::solver(" = TIME STEP DECREASED TO NEW VALUE                                                 %8.6f =\n", _nTimeShift);
			}
			eslog::solver(" = INCREASE FACTOR                                                                         %5.2f =\n", _nTimeShift / step::time.shift);
			changeConstants = true;
		}
		eslog::checkpointln("PHYSICS SOLVER: TIME STEP COMPUTED");
	}

	if (step::time.shift - step::time.precision < _nTimeShift) {
		Z->sum( _newmarkConsts[0], dU, -_newmarkConsts[2], V);
		Z->add(-_newmarkConsts[3], W);
		V->add( _newmarkConsts[6], W);
		V->add( _newmarkConsts[7], Z);
		std::swap(W, Z);

		U->fillData(system->solver()->x);
		system->processSolution();
	} else {
		system->solver()->x->fillData(U);
		system->solutionChanged();
		step::time.current = step::time.previous;
		--step::step.substep;
	}

	if (changeConstants) {
		updateConstants();
	}

	eslog::solver(" = ========================================================================================= =\n");
	eslog::solver(" = ================================================================= run time %12.3f s =\n\n", eslog::duration());
}


