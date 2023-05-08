

#include "analysis.h"
#include "structuralmechanics.steadystate.nonlinear.h"

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

StructuralMechanicsSteadyStateNonLinear::StructuralMechanicsSteadyStateNonLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, solver(configuration.nonlinear_solver), scheme{}, system{}
{

}

StructuralMechanicsSteadyStateNonLinear::~StructuralMechanicsSteadyStateNonLinear()
{
	if (system) {
		delete system;
	}
}

void StructuralMechanicsSteadyStateNonLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                        PLASTICITY STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                            STRUCTURAL MECHANICS == \n");
	eslog::info(" ============================================================================================= \n");

	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void StructuralMechanicsSteadyStateNonLinear::run(step::Step &step)
{
	initSystem(system, this);
	solver.init(system);
	scheme.init(system);
	assembler.connect(scheme.K, nullptr, nullptr, scheme.f, solver.R, scheme.x, scheme.dirichlet);
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

	solver.run(step, time, assembler, scheme, system);
	info::mesh->output->updateSolution(step, time);
}
