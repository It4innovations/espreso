
#include "looper.h"

#include "analysis/acoustic.real.linear.h"
#include "analysis/acoustic.complex.linear.h"
#include "analysis/heat.steadystate.linear.h"
#include "analysis/heat.steadystate.nonlinear.h"
#include "analysis/elasticity.steadystate.linear.h"

#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "output/output.h"

using namespace espreso;

void Looper::run()
{
	eslog::startln("ESPRESO: SIMULATION STARTED", "SIMULATION");

//	for (auto range = info::ecf->ranges.begin(); range != info::ecf->ranges.end(); ++range) {
//
//	}

	step::Step step;
	Analysis *analysis;

	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::ACOUSTICS_2D:
		switch (info::ecf->acoustics_2d.load_steps_settings.at(1).system) {
		case AcousticLoadStepConfiguration::SYSTEM::REAL: analysis = new AcousticRealLinear(info::ecf->acoustics_2d, info::ecf->acoustics_2d.load_steps_settings.at(1)); break;
		case AcousticLoadStepConfiguration::SYSTEM::COMPLEX: analysis = new AcousticComplexLinear(info::ecf->acoustics_2d, info::ecf->acoustics_2d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::ACOUSTICS_3D:
		switch (info::ecf->acoustics_3d.load_steps_settings.at(1).system) {
			case AcousticLoadStepConfiguration::SYSTEM::REAL: analysis = new AcousticRealLinear(info::ecf->acoustics_3d, info::ecf->acoustics_3d.load_steps_settings.at(1)); break;
			case AcousticLoadStepConfiguration::SYSTEM::COMPLEX: analysis = new AcousticComplexLinear(info::ecf->acoustics_3d, info::ecf->acoustics_3d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		switch (info::ecf->heat_transfer_2d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR: analysis = new HeatSteadyStateLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
		case LoadStepSolverConfiguration::MODE::NONLINEAR: analysis = new HeatSteadyStateNonLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		switch (info::ecf->heat_transfer_3d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR: analysis = new HeatSteadyStateLinear(info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(1)); break;
		case LoadStepSolverConfiguration::MODE::NONLINEAR: analysis = new HeatSteadyStateNonLinear(info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
		switch (info::ecf->structural_mechanics_2d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR: analysis = new ElasticitySteadyStateLinear(info::ecf->structural_mechanics_2d, info::ecf->structural_mechanics_2d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		switch (info::ecf->structural_mechanics_3d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR: analysis = new ElasticitySteadyStateLinear(info::ecf->structural_mechanics_3d, info::ecf->structural_mechanics_3d.load_steps_settings.at(1)); break;
		}
		break;
	}

	if (info::ecf->dryrun) {
		analysis->analyze();
		analysis->dryrun();
		// analysis->dryrun();
		// analysis->dryrun();
		eslog::checkpointln("SIMULATION: ASSEMBLER TESTED");
	} else {
		analysis->analyze();
		eslog::checkpointln("SIMULATION: PHYSICS ANALYSED");
		analysis->run(step);
	}

	delete analysis;
	eslog::endln("SIMULATION: DATA CLEARED");
}
