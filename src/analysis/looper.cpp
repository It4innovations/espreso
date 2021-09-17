
#include "looper.h"

#include "analysis/acoustic.real.linear.h"
#include "analysis/acoustic.complex.linear.h"
#include "analysis/heat.steadystate.linear.h"
#include "analysis/heat.steadystate.nonlinear.h"

#include "basis/expression/expression.h"
#include "basis/expression/variable.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "output/output.h"

using namespace espreso;

void Looper::run()
{
	Variable::gather(info::mesh->boundaryRegions.size());

//	for (auto range = info::ecf->ranges.begin(); range != info::ecf->ranges.end(); ++range) {
//
//	}

	step::Step step;
	Analysis *analysis;

	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::ACOUSTICS_2D:
		switch (info::ecf->acoustics_2d.load_steps_settings.at(1).system) {
		case AcousticLoadStepConfiguration::SYSTEM::REAL: analysis = new AX_AcousticRealLinear(info::ecf->acoustics_2d, info::ecf->acoustics_2d.load_steps_settings.at(1)); break;
		case AcousticLoadStepConfiguration::SYSTEM::COMPLEX: analysis = new AX_AcousticComplexLinear(info::ecf->acoustics_2d, info::ecf->acoustics_2d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		switch (info::ecf->heat_transfer_2d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR : analysis = new AX_HeatSteadyStateLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
		case LoadStepSolverConfiguration::MODE::NONLINEAR : analysis = new AX_HeatSteadyStateNonLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
		}
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		switch (info::ecf->heat_transfer_3d.load_steps_settings.at(1).mode) {
		case LoadStepSolverConfiguration::MODE::LINEAR : analysis = new AX_HeatSteadyStateNonLinear(info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(1)); break;
		case LoadStepSolverConfiguration::MODE::NONLINEAR : analysis = new AX_HeatSteadyStateNonLinear(info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(1)); break;
		}
		break;
	}

	analysis->init();
	analysis->run(step);

	Variable::clear();
}
