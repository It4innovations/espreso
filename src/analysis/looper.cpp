
#include "looper.h"

#include "analysis/acoustic.real.linear.h"
#include "analysis/heat.steadystate.linear.h"

#include "basis/expression/expression.h"
#include "basis/expression/variable.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "output/output.h"

using namespace espreso;

void Looper::run()
{
	Variable::gather(info::mesh->boundaryRegions.size());

//	for (auto range = info::ecf->ranges.begin(); range != info::ecf->ranges.end(); ++range) {
//
//	}

	Analysis *analysis;

	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::ACOUSTIC_2D:      analysis = new AX_AcousticRealLinear   (info::ecf->acoustic_2d     , info::ecf->acoustic_2d.load_steps_settings.at(1)); break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D: analysis = new AX_HeatSteadyStateLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D: analysis = new AX_HeatSteadyStateLinear(info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(1)); break;
	}

	analysis->init();
	analysis->run();
}
