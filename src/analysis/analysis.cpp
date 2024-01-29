
#include "analysis.h"
#include "physics/acoustic.real.linear.h"
#include "physics/acoustic.complex.linear.h"
#include "physics/heat.steadystate.linear.h"
#include "physics/heat.steadystate.nonlinear.h"
#include "physics/heat.transient.linear.h"
#include "physics/structuralmechanics.steadystate.linear.h"
#include "physics/structuralmechanics.steadystate.nonlinear.h"
#include "physics/structuralmechanics.harmonic.real.linear.h"

#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "output/output.h"

using namespace espreso;

void Analysis::run()
{
	eslog::startln("ESPRESO: SIMULATION STARTED", "SIMULATION");

//	for (auto range = info::ecf->ranges.begin(); range != info::ecf->ranges.end(); ++range) {
//
//	}

	step::Step step;
	Physics *physics;

	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::ACOUSTICS_2D:
		switch (info::ecf->acoustics_2d.load_steps_settings.at(1).system) {
		case AcousticLoadStepConfiguration::SYSTEM::REAL: physics = new AcousticRealLinear(info::ecf->acoustics_2d, info::ecf->acoustics_2d.load_steps_settings.at(1)); break;
		case AcousticLoadStepConfiguration::SYSTEM::COMPLEX: physics = new AcousticComplexLinear(info::ecf->acoustics_2d, info::ecf->acoustics_2d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::ACOUSTICS_3D:
		switch (info::ecf->acoustics_3d.load_steps_settings.at(1).system) {
			case AcousticLoadStepConfiguration::SYSTEM::REAL: physics = new AcousticRealLinear(info::ecf->acoustics_3d, info::ecf->acoustics_3d.load_steps_settings.at(1)); break;
			case AcousticLoadStepConfiguration::SYSTEM::COMPLEX: physics = new AcousticComplexLinear(info::ecf->acoustics_3d, info::ecf->acoustics_3d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		switch (info::ecf->heat_transfer_2d.load_steps_settings.at(1).type) {
		case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
			switch (info::ecf->heat_transfer_2d.load_steps_settings.at(1).mode) {
			case LoadStepSolverConfiguration::MODE::LINEAR: physics = new HeatSteadyStateLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
			case LoadStepSolverConfiguration::MODE::NONLINEAR: physics = new HeatSteadyStateNonLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
			} break;
		case LoadStepSolverConfiguration::TYPE::TRANSIENT:
			switch (info::ecf->heat_transfer_2d.load_steps_settings.at(1).mode) {
			case LoadStepSolverConfiguration::MODE::LINEAR: physics = new HeatTransientLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
//			case LoadStepSolverConfiguration::MODE::NONLINEAR: physics = new HeatSteadyStateNonLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
			} break;
		}

		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		switch (info::ecf->heat_transfer_3d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR: physics = new HeatSteadyStateLinear(info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(1)); break;
		case LoadStepSolverConfiguration::MODE::NONLINEAR: physics = new HeatSteadyStateNonLinear(info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
		switch (info::ecf->structural_mechanics_2d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR: physics = new StructuralMechanicsSteadyStateLinear(info::ecf->structural_mechanics_2d, info::ecf->structural_mechanics_2d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		switch (info::ecf->structural_mechanics_3d.load_steps_settings.at(1).type) {
		case LoadStepSolverConfiguration::TYPE::HARMONIC:
			switch (info::ecf->structural_mechanics_3d.load_steps_settings.at(1).mode) {
			case LoadStepSolverConfiguration::MODE::LINEAR: physics = new StructuralMechanicsHarmonicRealLinear(info::ecf->structural_mechanics_3d, info::ecf->structural_mechanics_3d.load_steps_settings.at(1)); break;
			case LoadStepSolverConfiguration::MODE::NONLINEAR:  break;
			} break;
		case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
			switch (info::ecf->structural_mechanics_3d.load_steps_settings.at(1).mode) {
			case LoadStepSolverConfiguration::MODE::LINEAR: physics = new StructuralMechanicsSteadyStateLinear(info::ecf->structural_mechanics_3d, info::ecf->structural_mechanics_3d.load_steps_settings.at(1)); break;
			case LoadStepSolverConfiguration::MODE::NONLINEAR: physics = new StructuralMechanicsSteadyStateNonLinear(info::ecf->structural_mechanics_3d, info::ecf->structural_mechanics_3d.load_steps_settings.at(1)); break;
			} break;
		case LoadStepSolverConfiguration::TYPE::TRANSIENT:
			break;
		} break;
		break;
	}

	physics->analyze(step);
	eslog::checkpointln("SIMULATION: PHYSICS ANALYSED");
	step.loadstep = 0;
	step.loadsteps = 1;
	physics->run(step);

	delete physics;
	eslog::endln("SIMULATION: DATA CLEARED");
}
