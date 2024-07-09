
#include "physics.h"
#include "heat.steadystate.linear.h"

#include "analysis/builder/uniformbuilder.direct.h"
#include "analysis/builder/uniformbuilder.feti.h"
#include "analysis/linearsystem/mklpdsssolver.h"
#include "analysis/linearsystem/fetisolver.h"
#include "analysis/linearsystem/empty.h"
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
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, f{}, x{}, dirichlet{}, builder{}, solver{}
{

}

HeatSteadyStateLinear::~HeatSteadyStateLinear()
{
    if (K) { delete K; }
    if (f) { delete f; }
    if (x) { delete x; }
    if (dirichlet) { delete dirichlet; }
    if (builder) { delete builder; }
    if (solver) { delete solver; }
}

void HeatSteadyStateLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == PHYSICS                                                                   HEAT TRANSFER == \n");
    eslog::info(" == ANALYSIS                                                                   STEADY STATE == \n");
    eslog::info(" == MODE                                                                             LINEAR == \n");
    eslog::info(" ============================================================================================= \n");

    step.type = step::TYPE::TIME;
    assembler.analyze();
    info::mesh->output->updateMonitors(step);

    switch (configuration.solver) {
    case LoadStepSolverConfiguration::SOLVER::FETI:
        builder = new UniformBuilderFETI<double>(configuration);
        solver = new FETILinearSystemSolver<double>(settings, configuration);
        break;
    case LoadStepSolverConfiguration::SOLVER::HYPRE:   break;
    case LoadStepSolverConfiguration::SOLVER::MKLPDSS:
        builder = new UniformBuilderDirect<double>(configuration);
        solver = new MKLPDSSLinearSystemSolver<double>(configuration.mklpdss);
        break;
    case LoadStepSolverConfiguration::SOLVER::PARDISO: break;
    case LoadStepSolverConfiguration::SOLVER::SUPERLU: break;
    case LoadStepSolverConfiguration::SOLVER::WSMP:    break;
    case LoadStepSolverConfiguration::SOLVER::NONE:
        builder = new UniformBuilderDirect<double>(configuration);
        solver = new EmptySystemSolver<double>();
    }

    builder->fillMatrix(solver->A);
    builder->fillVector(solver->b);
    builder->fillVector(solver->x);
    builder->fillDirichlet(solver->dirichlet);

    K = solver->A->copyPattern();
    f = solver->b->copyPattern();
    x = solver->x->copyPattern();
    dirichlet = solver->dirichlet->copyPattern();

    builder->fillMatrixMap(K);
    builder->fillVectorMap(f);
    builder->fillDirichletMap(dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
}

void HeatSteadyStateLinear::run(step::Step &step)
{
    time.shift = configuration.duration_time;
    time.start = 0;
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
    solver->solve(step);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");
    double solution = eslog::time();

    x->copy(solver->x);
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

