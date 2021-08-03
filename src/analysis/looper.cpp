
#include "looper.h"

#include "analysis/acoustic.real.linear.h"
#include "analysis/heat.steadystate.linear.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "output/output.h"

#include <cstdio>

using namespace espreso;

void Looper::run()
{
	Analysis *analysis;

	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::ACOUSTIC_2D:      analysis = new AX_AcousticRealLinear   (info::ecf->acoustic_2d     , info::ecf->acoustic_2d.load_steps_settings.at(1)); break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D: analysis = new AX_HeatSteadyStateLinear(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1)); break;
	}

	analysis->init();
	analysis->run();
}
