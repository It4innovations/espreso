
#include "analysis.h"
#include "heat.steadystate.linear.h"

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

HeatSteadyStateLinear::HeatSteadyStateLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, f{}, x{}, dirichlet{}, system{}
{

}

HeatSteadyStateLinear::~HeatSteadyStateLinear()
{
	if (system) { delete system; }
	if (K) { delete K; }
	if (f) { delete f; }
	if (x) { delete x; }
	if (dirichlet) { delete dirichlet; }
}

void HeatSteadyStateLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                            LINEAR STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                                   HEAT TRANSFER == \n");
	eslog::info(" ============================================================================================= \n");

	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void HeatSteadyStateLinear::run(step::Step &step)
{
	initSystem(system, this);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");

	system->setMapping(K = system->assembler.A->copyPattern());
	system->setMapping(f = system->assembler.b->copyPattern());
	system->setMapping(x = system->assembler.x->copyPattern());
	system->setDirichletMapping(dirichlet = system->assembler.dirichlet->copyPattern());
	assembler.connect(K, nullptr, f, nullptr, dirichlet);

	time.shift = configuration.duration_time;
	time.start = 0;
	time.current = configuration.duration_time;
	time.final = configuration.duration_time;

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
	assembler.evaluate(time, K, nullptr, f, nullptr, dirichlet);
	eslog::checkpointln("SIMULATION: PHYSICS ASSEMBLED");
	storeSystem(step);

	system->solver.A->copy(K);
	system->solver.b->copy(f);
	system->solver.dirichlet->copy(dirichlet);

	eslog::info("       = ----------------------------------------------------------------------------- = \n");
	eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

	system->update(step);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM UPDATED");
	system->solve(step);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");
	double solution = eslog::time();

	x->copy(system->solver.x);
	storeSolution(step);
	assembler.updateSolution(x);
	info::mesh->output->updateSolution(step, time);
	eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
	eslog::info("       = ----------------------------------------------------------------------------- = \n");

	eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
	eslog::checkpointln("SIMULATION: SOLUTION PROCESSED");
}

void HeatSteadyStateLinear::storeSystem(step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{K, f}\n");
		K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
		f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
		dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet").c_str());
	}
}

void HeatSteadyStateLinear::storeSolution(step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{x}\n");
		x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
	}
}

