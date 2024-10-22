
#include "structuralmechanics.steadystate.linear.h"

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
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, f{}, x{}, forces{}, dirichlet{}, pattern{}, solver{}
{

}

StructuralMechanicsSteadyStateLinear::~StructuralMechanicsSteadyStateLinear()
{
    if (K) { delete K; }
    if (f) { delete f; }
    if (x) { delete x; }
    if (forces) { delete forces; }
    if (dirichlet) { delete dirichlet; }
    if (pattern) { delete pattern; }
    if (solver) { delete solver; }
}

bool StructuralMechanicsSteadyStateLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == ANALYSIS                                                                   STEADY STATE == \n");
    eslog::info(" == PHYSICS                                                            STRUCTURAL MECHANICS == \n");
    eslog::info(" == MODE                                                                             LINEAR == \n");
    eslog::info(" ============================================================================================= \n");

    step.type = step::TYPE::TIME;
    if (!assembler.analyze(step)) {
        return false;
    }
    info::mesh->output->updateMonitors(step);

    solver = setSolver<double>(settings, configuration);
    pattern = solver->getPattern(configuration, 1);

    pattern->set(solver->A);
    pattern->set(solver->b);
    pattern->set(solver->x);
    pattern->set(solver->dirichlet);

    K = solver->A->copyPattern();
    f = solver->b->copyPattern();
    x = solver->x->copyPattern();
    dirichlet = solver->dirichlet->copyPattern();
    forces = solver->x->copyPattern();

    pattern->map(K);
    pattern->map(f);
    pattern->map(dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
    return true;
}

bool StructuralMechanicsSteadyStateLinear::run(step::Step &step, Physics *prev)
{
    step.substep = 0;
    step.substeps = 1;
    time.shift = configuration.duration_time;
    time.start = 0;
    if (prev) {
        bool correct = false;
        if (!correct) {
            eslog::globalerror("Incompatible load steps.\n");
        }
    }
    time.current = configuration.duration_time;
    time.final = configuration.duration_time;

    assembler.connect(K, nullptr, f, nullptr, dirichlet);

    if (MPITools::node->rank == 0) {
        info::system::memory::physics = info::system::memoryAvail();
    }
    eslog::info("  PHYSICAL SOLVER MEMORY FOOTPRINT [GB] %53.2f  \n", (info::system::memory::mesh - info::system::memory::physics) / 1024. / 1024.);
    eslog::info(" ============================================================================================= \n");

    eslog::info("\n ============================================================================================= \n");
    eslog::info(" = RUN THE SOLVER                                                DURATION TIME: %10.4f s = \n", configuration.duration_time);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    solver->set(step);
    eslog::info(" ============================================================================================= \n\n");
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM SET");

    eslog::info(" ============================================================================================= \n");
    eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step.loadstep + 1, time.current);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    double start = eslog::time();
    assembler.evaluate(step, time, K, nullptr, f, nullptr, dirichlet);
    eslog::checkpointln("SIMULATION: PHYSICS ASSEMBLED");
    storeSystem(step);

    solver->A->copy(K);
    solver->A->updated = true;
    solver->b->copy(f);
    solver->dirichlet->copy(dirichlet);
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
    eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

    solver->update(step);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM UPDATED");
    bool solved = solver->solve(step);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");

    double solution = eslog::time();
    x->copy(solver->x);
    if (StructuralMechanics::Results::reactionForce) {
        K->apply(1., x, 0., forces);
        forces->storeTo(StructuralMechanics::Results::reactionForce->data);
    }
    storeSolution(step);
    assembler.updateSolution(x);
    info::mesh->output->updateSolution(step, time);

    eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
    eslog::info("       = ----------------------------------------------------------------------------- = \n");

    eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
    eslog::checkpointln("SIMULATION: SOLUTION PROCESSED");
    return solved;
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
