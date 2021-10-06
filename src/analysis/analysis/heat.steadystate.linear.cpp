
#include "analysis.h"
#include "heat.steadystate.linear.h"

#include "analysis/linearsystem/linearsystem.hpp"
#include "basis/expression/variable.h"
#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"

using namespace espreso;

AX_HeatSteadyStateLinear::AX_HeatSteadyStateLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, scheme{}, system{}
{

}

AX_HeatSteadyStateLinear::~AX_HeatSteadyStateLinear()
{
	if (system) {
		delete system;
	}
}

void AX_HeatSteadyStateLinear::init()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                            LINEAR STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                                   HEAT TRANSFER == \n");
	eslog::info(" ============================================================================================= \n");

	initSystem(system, this);
	scheme.init(system);
	assembler.init(scheme);

	Variable::list.global.insert(std::make_pair("TIME", nullptr));

	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void AX_HeatSteadyStateLinear::run(step::Step &step)
{
	step::Time time;
	scheme.setTime(time, configuration.duration_time);
	Variable::list.global["TIME"] = new TimeVariable(time);

	eslog::info("\n ============================================================================================= \n");
	eslog::info(" = RUN THE SOLVER                                                DURATION TIME: %10.4f s = \n", configuration.duration_time);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	system->info();
	system->set(step);
	eslog::info(" ============================================================================================= \n\n");

	eslog::info(" ============================================================================================= \n");
	eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step::step.loadstep + 1, time.current);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	double start = eslog::time();
	assembler.evaluate();
	scheme.composeSystem(step, system);
	eslog::info("       = ----------------------------------------------------------------------------- = \n");
	eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

	system->update(step);
	system->solve(step);

	double solution = eslog::time();
	scheme.extractSolution(step, system);
	assembler.updateSolution();
	info::mesh->output->updateSolution(step, time);
	eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
	eslog::info("       = ----------------------------------------------------------------------------- = \n");

	eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
}


