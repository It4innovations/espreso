
#include "looper.h"

#include "analysis/heat.steadystate.linear.h"

#include "esinfo/ecfinfo.h"

#include <cstdio>

using namespace espreso;

void Looper::run()
{
	printf("looper\n");

	AX_HeatSteadyStateLinear analysis(info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(1));
}
