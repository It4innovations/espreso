
#include "advectiondiffusionfactory.h"

#include "../../configuration/physics/advectiondiffusion2d.h"
#include "../../configuration/physics/advectiondiffusion3d.h"

#include "../../assembler/physics/advectiondiffusion2d.h"
#include "../../assembler/physics/advectiondiffusion3d.h"
#include "../../assembler/physics/laplacesteklovpoincare3d.h"

#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physicssolver/timestep/linear.h"
#include "../../assembler/physicssolver/timestep/newtonrhapson.h"
#include "../../assembler/physicssolver/loadstep/steadystate.h"
#include "../../assembler/physicssolver/loadstep/pseudotimestepping.h"
#include "../../assembler/physicssolver/loadstep/transientfirstorderimplicit.h"

#include "../../assembler/instance.h"
#include "../../mesh/structures/mesh.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

AdvectionDiffusionFactory::AdvectionDiffusionFactory(const AdvectionDiffusion2DConfiguration &configuration, Mesh *mesh)
: _configuration(configuration), _bem(false)
{
	_instances.push_back(new Instance(mesh->parts(), mesh->neighbours()));
	_physics.push_back(new AdvectionDiffusion2D(mesh, _instances.front(), configuration));
}

AdvectionDiffusionFactory::AdvectionDiffusionFactory(const AdvectionDiffusion3DConfiguration &configuration, Mesh *mesh)
: _configuration(configuration), _bem(false)
{
	_instances.push_back(new Instance(mesh->parts(), mesh->neighbours()));

	switch (configuration.discretization) {
	case DISCRETIZATION::FEM:
		_physics.push_back(new AdvectionDiffusion3D(mesh, _instances.front(), configuration));
		break;
	case DISCRETIZATION::BEM:
		_bem = true;
		_physics.push_back(new LaplaceSteklovPoincare3D(mesh, _instances.front(), configuration));
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown DISCRETIZATION for AdvectionDiffusion3D";
	}
}

size_t AdvectionDiffusionFactory::loadSteps() const
{
	return _configuration.physics_solver.load_steps;
}

LoadStepSolver* AdvectionDiffusionFactory::getLoadStepSolver(size_t step, Mesh *mesh, output::Store *store)
{
	const LoadStepSettings<AdvectionDiffusionNonLinearConvergence> &settings = getLoadStepsSettings(step, _configuration.physics_solver.load_steps_settings);

	_linearSolvers.push_back(getLinearSolver(settings, _instances.front()));
	_assemblers.push_back(new Assembler(*_instances.front(), *_physics.front(), *mesh, *store, *_linearSolvers.back()));

	switch (settings.mode) {
	case LoadStepSettingsBase::MODE::LINEAR:
		_timeStepSolvers.push_back(new LinearTimeStep(*_assemblers.back()));
		break;
	case LoadStepSettingsBase::MODE::NONLINEAR:
		if (_bem) {
			ESINFO(GLOBAL_ERROR) << "BEM discretization support only LINEAR STEADY STATE physics solver.";
		}
		switch (settings.nonlinear_solver.method) {
		case NonLinearSolverBase::METHOD::NEWTON_RHAPSON:
		case NonLinearSolverBase::METHOD::MODIFIED_NEWTON_RHAPSON:
			_timeStepSolvers.push_back(new NewtonRhapson(*_assemblers.back(), settings.nonlinear_solver));
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented NONLINEAR SOLVER METHOD for LOAD STEP=" << step;
		}
		break;

	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented LOAD STEP solver MODE for LOAD STEP=" << step;
	}

	switch (settings.type) {
	case LoadStepSettingsBase::TYPE::STEADY_STATE:
		if (settings.mode == LoadStepSettingsBase::MODE::NONLINEAR) {
			_loadStepSolvers.push_back(new PseudoTimeStepping(*_timeStepSolvers.back(), settings.nonlinear_solver, settings.duration_time));
		} else {
			_loadStepSolvers.push_back(new SteadyStateSolver(*_timeStepSolvers.back(), settings.duration_time));
		}
		break;
	case LoadStepSettingsBase::TYPE::TRANSIENT:
		_loadStepSolvers.push_back(new TransientFirstOrderImplicit(*_timeStepSolvers.back(), settings.transient_solver, settings.duration_time));
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented LOAD STEP solver TYPE for LOAD STEP=" << step;
	}

	return _loadStepSolvers.back();
}


