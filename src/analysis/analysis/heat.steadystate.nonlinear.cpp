
#include "analysis.h"
#include "heat.steadystate.nonlinear.h"

#include "analysis/linearsystem/linearsystem.hpp"
#include "basis/expression/variable.h"
#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"

using namespace espreso;

HeatSteadyStateNonLinear::HeatSteadyStateNonLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, solver(configuration.nonlinear_solver), scheme{}, system{}

{

}

HeatSteadyStateNonLinear::~HeatSteadyStateNonLinear()
{
	if (system) {
		delete system;
	}
}

void HeatSteadyStateNonLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                        NON-LINEAR STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                                   HEAT TRANSFER == \n");
	eslog::info(" ============================================================================================= \n");

	Variable::list.global.insert(std::make_pair("TIME", new TimeVariable(time)));
	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void HeatSteadyStateNonLinear::run(step::Step &step)
{
	initSystem(system, this);
	solver.init(system);
	scheme.init(system);
	assembler.connect(scheme);
	scheme.setTime(time, configuration.duration_time);

	eslog::info("\n ============================================================================================= \n");
	eslog::info(" = RUN THE SOLVER                                                DURATION TIME: %10.4f s = \n", configuration.duration_time);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	system->set(step);
	eslog::info(" ============================================================================================= \n\n");

	eslog::info(" ============================================================================================= \n");
	eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step::step.loadstep + 1, time.current);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");

	solver.run(step, time, assembler, scheme, system);

	info::mesh->output->updateSolution(step, time);
}



