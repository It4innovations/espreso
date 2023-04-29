
#include "analysis.h"
#include "heat.steadystate.nonlinear.h"

#include "analysis/linearsystem/linearsystem.hpp"
#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"

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
	if (MPITools::node->rank == 0) {
		info::system::memory::physics = info::system::memoryAvail();
	}
	eslog::info("  PHYSICAL SOLVER MEMORY FOOTPRINT [GB] %53.2f  \n", (info::system::memory::mesh - info::system::memory::physics) / 1024. / 1024.);
	eslog::info(" ============================================================================================= \n");

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



