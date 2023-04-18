
#include "analysis.h"
#include "elasticity.steadystate.linear.h"

#include "analysis/linearsystem/linearsystem.hpp"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

ElasticitySteadyStateLinear::ElasticitySteadyStateLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, scheme{}, system{}
{

}

ElasticitySteadyStateLinear::~ElasticitySteadyStateLinear()
{
	if (system) {
		delete system;
	}
}

void ElasticitySteadyStateLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                            LINEAR STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                            STRUCTURAL MECHANICS == \n");
	eslog::info(" ============================================================================================= \n");

	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void ElasticitySteadyStateLinear::run(step::Step &step)
{
	initSystem(system, this);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
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
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM SET");

	eslog::info(" ============================================================================================= \n");
	eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step::step.loadstep + 1, time.current);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	double start = eslog::time();
	assembler.evaluate(scheme);
	eslog::checkpointln("SIMULATION: PHYSICS ASSEMBLED");
	scheme.composeSystem(step, system);
	eslog::info("       = ----------------------------------------------------------------------------- = \n");
	eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

	system->update(step);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM UPDATED");
	system->solve(step);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");

	double solution = eslog::time();
	scheme.extractSolution(step, system);
	assembler.updateSolution(scheme);
	info::mesh->output->updateSolution(step, time);
	eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
	eslog::info("       = ----------------------------------------------------------------------------- = \n");

	eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
	eslog::checkpointln("SIMULATION: SOLUTION PROCESSED");
}
