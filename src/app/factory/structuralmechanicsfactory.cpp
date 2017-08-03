
#include "structuralmechanicsfactory.h"

#include "../../configuration/physics/structuralmechanics2d.h"
#include "../../configuration/physics/structuralmechanics3d.h"

#include "../../assembler/physics/structuralmechanics2d.h"
#include "../../assembler/physics/structuralmechanics3d.h"
#include "../../assembler/physics/lamesteklovpoincare3d.h"

#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physicssolver/timestep/linear.h"
#include "../../assembler/physicssolver/timestep/newtonrhapson.h"
#include "../../assembler/physicssolver/loadstep/steadystate.h"

#include "../../assembler/instance.h"
#include "../../mesh/structures/mesh.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

StructuralMechanicsFactory::StructuralMechanicsFactory(const StructuralMechanics2DConfiguration &configuration, Mesh *mesh)
: _configuration(configuration), _bem(false)
{
	_instances.push_back(new Instance(mesh->parts(), mesh->neighbours()));
	_physics.push_back(new StructuralMechanics2D(mesh, _instances.front(), configuration));
}

StructuralMechanicsFactory::StructuralMechanicsFactory(const StructuralMechanics3DConfiguration &configuration, Mesh *mesh)
: _configuration(configuration), _bem(false)
{
	_instances.push_back(new Instance(mesh->parts(), mesh->neighbours()));

	switch (configuration.discretization) {
	case DISCRETIZATION::FEM:
		_physics.push_back(new StructuralMechanics3D(mesh, _instances.front(), configuration));
		break;
	case DISCRETIZATION::BEM:
		_bem = true;
		_physics.push_back(new LameSteklovPoincare3D(mesh, _instances.front(), configuration));
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown DISCRETIZATION for StructuralMechanics3D";
	}
}

size_t StructuralMechanicsFactory::loadSteps() const
{
	return _configuration.physics_solver.load_steps;
}

LoadStepSolver* StructuralMechanicsFactory::getLoadStepSolver(size_t step, Mesh *mesh, output::Store *store)
{
	const LoadStepSettings<StructuralMechanicsNonLinearConvergence> &settings = getLoadStepsSettings(step, _configuration.physics_solver.load_steps_settings);

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
		_loadStepSolvers.push_back(new SteadyStateSolver(*_timeStepSolvers.back(), settings.duration_time));
		break;
	case LoadStepSettingsBase::TYPE::TRANSIENT:
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented LOAD STEP solver TYPE for LOAD STEP=" << step;
	}

	return _loadStepSolvers.back();
}



