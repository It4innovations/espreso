
#include "harmonicsolver.h"
#include "esinfo/stepinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"

#include "basis/containers/tarray.h"
#include "basis/evaluator/evaluator.h"

#include "physics/system/linearsystem.h"
#include "physics/system/builder/builder.h"
#include "physics/substepsolver/substepsolver.h"
#include "config/ecf/physics/physicssolver/loadstep.h"
#include "output/output.h"

#include <cmath>
#include <algorithm>

using namespace espreso;

HarmonicSolver::HarmonicSolver(LinearSystem *system, SubStepSolver *subStepSolver, HarmonicSolverConfiguration &configuration, double duration)
: LoadStepSolver(system, subStepSolver, duration), _configuration(configuration)
{
	std::vector<int> distribution = tarray<int>::distribute(info::mpi::isize, _configuration.num_samples);
	step::duplicate::instances = info::mpi::isize;
	step::duplicate::offset = distribution[info::mpi::irank];
	step::duplicate::size = distribution[info::mpi::irank + 1] - distribution[info::mpi::irank];
	step::duplicate::totalsize = _configuration.num_samples;

	step::type = step::TYPE::FREQUENCY;
	step::frequency::shift = (_configuration.max_frequency - _configuration.min_frequency) / _configuration.num_samples;
	step::frequency::start = _configuration.min_frequency;
	step::frequency::current = step::frequency::start + step::duplicate::offset * step::frequency::shift;
	step::frequency::final = _configuration.max_frequency;

	_fttRequestedFrequencies = info::ecf->output.frequency_to_time.requested_frequencies.data();
	_fttRequestedFrequenciesEnd = _fttRequestedFrequencies + info::ecf->output.frequency_to_time.requested_frequencies.size();

	std::sort(_fttRequestedFrequencies, _fttRequestedFrequenciesEnd);
	_fttRequestedFrequencies = std::lower_bound(_fttRequestedFrequencies, _fttRequestedFrequenciesEnd, step::frequency::current);

	updateDamping();
}

void HarmonicSolver::init(LoadStepSolver *previous)
{

}

void HarmonicSolver::updateDamping()
{
	switch (_configuration.damping.rayleigh.type) {
	case RayleighDampingConfiguration::TYPE::NONE:
		_system->builder->rayleighDamping = false;
		_system->builder->stiffnessDamping = 0;
		_system->builder->massDamping = 0;
		break;
	case RayleighDampingConfiguration::TYPE::DIRECT:
		_system->builder->rayleighDamping = true;
		_system->builder->stiffnessDamping = _configuration.damping.rayleigh.direct_damping.stiffness.evaluator->eval({});
		_system->builder->massDamping = _configuration.damping.rayleigh.direct_damping.mass.evaluator->eval({});
		break;
	case RayleighDampingConfiguration::TYPE::DAMPING_RATIO: {
		_system->builder->rayleighDamping = true;
		double ratio = _configuration.damping.rayleigh.ratio_damping.ratio.evaluator->eval({});
		double frequency = _configuration.damping.rayleigh.ratio_damping.frequency.evaluator->eval({});
		_system->builder->stiffnessDamping = 2 * ratio * 2 * M_PI * frequency;
		_system->builder->massDamping = 2 * ratio / (2 * M_PI * frequency);
	} break;
	}

	_system->builder->coriolisDamping = _configuration.damping.coriolis_effect.coriolis_damping;
	_system->builder->rotationAxis.x = _configuration.damping.coriolis_effect.rotation_axis.x;
	_system->builder->rotationAxis.y = _configuration.damping.coriolis_effect.rotation_axis.y;
	_system->builder->rotationAxis.z = _configuration.damping.coriolis_effect.rotation_axis.z;
	_system->builder->rotationAxis.normalize();
}

void HarmonicSolver::store()
{
	switch (info::ecf->output.frequency_to_time.results_store_frequency) {
	case HarmonicOuputConfiguration::STORE_FREQUENCY::NEVER:
		_system->processSolution(); break;
	case HarmonicOuputConfiguration::STORE_FREQUENCY::EVERY_FREQUENCY:
		ftt(); break;
	case HarmonicOuputConfiguration::STORE_FREQUENCY::EVERY_NTH_FREQUENCY:
		if (step::substep % info::ecf->output.frequency_to_time.results_nth_stepping == 0) {
			ftt();
		} else {
			_system->processSolution();
		}
		break;
	case HarmonicOuputConfiguration::STORE_FREQUENCY::SPECIFIC_FREQUENCIES:
		while (_fttRequestedFrequencies != _fttRequestedFrequenciesEnd && *_fttRequestedFrequencies < step::frequency::current) {
			++_fttRequestedFrequencies;
		}
		if (_fttRequestedFrequencies != _fttRequestedFrequenciesEnd && *_fttRequestedFrequencies == step::frequency::current) {
			ftt();
		} else {
			_system->processSolution();
		}
		break;
	}
}

void HarmonicSolver::ftt()
{
	_system->processSolution();
	step::type = step::TYPE::FTT;
	step::ftt::steps = info::ecf->output.frequency_to_time.samples;
	step::ftt::period = 1 / step::frequency::current;
	for (int i = 0; i < step::ftt::steps; i++) {
		step::ftt::step = i;
		step::ftt::time = step::ftt::period * ((double)i / step::ftt::steps);
		_system->processSolution();
	}
	step::type = step::TYPE::FREQUENCY;
}

void HarmonicSolver::updateStructuralMatrices()
{
	Builder::Request matrices = Builder::Request::K | Builder::Request::M | Builder::Request::RBCf;
	if (
			(_system->builder->rayleighDamping && _system->builder->spinSoftening) ||
			(_system->builder->coriolisDamping)) {

		matrices |= Builder::Request::C;
	}

	_system->builder->matrices &= matrices;
	_system->assemble();
}

void HarmonicSolver::runNextSubstep()
{
	switch (_configuration.frequency_interval_type) {
	case HarmonicSolverConfiguration::INTERVAL_TYPE::LINEAR:
		step::frequency::current += step::frequency::shift;
		if (step::frequency::current + step::frequency::precision >= step::frequency::final) {
			step::frequency::current = step::frequency::final;
		}
		break;
	default:
		eslog::internalFailure("not implemented interval type.\n");
	}
	step::frequency::angular = 2 * M_PI * step::frequency::current;
	_system->nextSubstep();

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

	_system->builder->internalForceReduction = 1;
	
	switch (_configuration.aft.type) {
	case AlternatingFrequencyTime::TYPE::USER:
		_system->builder->AFTSamples = _configuration.aft.time_samples;
		break;
	default:
		eslog::error("AFT type not implemented!\n");
	}

	eslog::solver(" = ==================================== HARMONIC SOLVER ==================================== =\n");
	eslog::solver(" =  LOAD STEP %2d, SUBSTEP %4d, FREQ %10.4f, FREQ STEP %10.4f, FINAL FREQ %10.4f =\n", step::loadstep + 1, step::substep + 1, step::frequency::current, step::frequency::shift, step::frequency::final);
	eslog::solver(" = ----------------------------------------------------------------------------------------- =\n");

	_subStepSolver->solve(*this);
	store();

	eslog::solver(" = ========================================================================================= =\n");
	eslog::solver(" = ================================================================= run time %12.3f s =\n\n", eslog::duration());
}

