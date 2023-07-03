
#include <analysis/physics/physics.h>
#include "structuralmechanics.steadystate.linear.h"

#include "analysis/linearsystem/feti/fetisystem.h"
#include "analysis/linearsystem/direct/mklpdsssystem.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

StructuralMechanicsSteadyStateLinear::StructuralMechanicsSteadyStateLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, f{}, x{}, dirichlet{}, system{}
{

}

StructuralMechanicsSteadyStateLinear::~StructuralMechanicsSteadyStateLinear()
{
	if (system) { delete system; }
	if (K) { delete K; }
	if (f) { delete f; }
	if (x) { delete x; }
	if (dirichlet) { delete dirichlet; }
}

void StructuralMechanicsSteadyStateLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                            LINEAR STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                            STRUCTURAL MECHANICS == \n");
	eslog::info(" ============================================================================================= \n");

	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void StructuralMechanicsSteadyStateLinear::run(step::Step &step)
{
	switch (configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new FETISystem<StructuralMechanicsSteadyStateLinear>(this); break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new MKLPDSSSystem<StructuralMechanicsSteadyStateLinear>(this); break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    break;
	}
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");

	system->setMapping(K = system->assembler.A->copyPattern());
	system->setMapping(f = system->assembler.b->copyPattern());
	system->setMapping(x = system->assembler.x->copyPattern());
	system->setDirichletMapping(dirichlet = system->assembler.dirichlet->copyPattern());
	assembler.connect(K, nullptr, nullptr, f, nullptr, dirichlet);

	step.substep = 0;
	step.substeps = 1;
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
	assembler.evaluate(step, time, K, nullptr, nullptr, f, nullptr, dirichlet);
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

void StructuralMechanicsSteadyStateLinear::storeSystem(step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{K, f}\n");
		K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
		f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
		dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet").c_str());
	}
}

void StructuralMechanicsSteadyStateLinear::storeSolution(step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{x}\n");
		x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
	}
}
