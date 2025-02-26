
#include "physics.h"
#include "heat.steadystate.nonlinear.h"

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

HeatSteadyStateNonLinear::HeatSteadyStateNonLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{settings, configuration}, K{}, U{}, R{}, f{}, x{}, dirichlet{}, pattern{}, solver{}

{

}

HeatSteadyStateNonLinear::~HeatSteadyStateNonLinear()
{
    if (K) { delete K; }
    if (U) { delete U; }
    if (R) { delete R; }
    if (f) { delete f; }
    if (x) { delete x; }
    if (dirichlet) { delete dirichlet; }
    if (pattern) { delete pattern; }
    if (solver) { delete solver; }
}

bool HeatSteadyStateNonLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == ANALYSIS                                                                   STEADY STATE == \n");
    eslog::info(" == PHYSICS                                                                   HEAT TRANSFER == \n");
    eslog::info(" == MODE                                                                         NON-LINEAR == \n");
    eslog::info(" ============================================================================================= \n");

    step.type = step::TYPE::TIME;
    if (!assembler.analyze(step)) {
        return false;
    }
    info::mesh->output->updateMonitors(step);

    solver = setSolver<double>(configuration);
    pattern = solver->getPattern(configuration, 1);

    pattern->set(solver);

    K = solver->A->copyPattern();
    R = solver->b->copyPattern();
    U = solver->b->copyPattern();
    f = solver->b->copyPattern();
    x = solver->x->copyPattern();
    dirichlet = solver->dirichlet->copyPattern();
    solver->assembledA = K;

    pattern->map(K);
    pattern->map(R);
    pattern->map(f);
    pattern->map(dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
    return true;
}

bool HeatSteadyStateNonLinear::run(step::Step &step, Physics *prev)
{
    time.start = 0;
    if (prev) {
        bool correct = false;
        if (!correct) {
            eslog::globalerror("Incompatible load steps.\n");
        }
    }

    time.shift = configuration.duration_time;
    time.current = time.start + time.shift;
    time.final   = time.start + configuration.duration_time;
    time.timeIntegrationConstantK = 1;

    assembler.connect(K, nullptr, f, R, dirichlet);
    R->set(0);

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

    eslog::info(" ============================================================================================= \n");
    eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step.loadstep + 1, time.current);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    eslog::info("      =================================================================================== \n");
    eslog::info("      ==  NEWTON RAPHSON CONVERGENCE CRITERIA                                          == \n");
    if (configuration.nonlinear_solver.check_first_residual) {
        eslog::info("      ==  - TEMPERATURE RESIDUAL                                                  TRUE == \n");
    } else {
        eslog::info("      ==  - TEMPERATURE RESIDUAL                                                 FALSE == \n");
    }
    if (configuration.nonlinear_solver.check_second_residual) {
        eslog::info("      ==  - HEAT RESIDUAL                                                         TRUE == \n");
    } else {
        eslog::info("      ==  - HEAT RESIDUAL                                                        FALSE == \n");
    }
    eslog::info("      =================================================================================== \n");

    eslog::info("      ==                                                                  INITIAL STEP ==     \n");

    double start = eslog::time();
    step.iteration = 0;
    assembler.evaluate(step, time, K, nullptr, f, nullptr, dirichlet);
    storeSystem(step);
    solver->A->copy(K);
    solver->A->updated = true;
    solver->b->copy(f);
    solver->dirichlet->copy(dirichlet);
    eslog::info("      == ----------------------------------------------------------------------------- == \n");
    eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

    solver->update(step);
    solver->solve(step);

    double solution = eslog::time();
    x->copy(solver->x);
    storeSolution(step);
    assembler.updateSolution(step, x);
    eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
    eslog::info("      == ----------------------------------------------------------------------------- == \n");

    bool converged = false;
    while (step.iteration++ < configuration.nonlinear_solver.max_iterations) {
        eslog::info("\n      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

        start = eslog::time();
        U->copy(solver->x);
        assembler.evaluate(step, time, K, nullptr, f, R, dirichlet);
        storeSystem(step);
        solver->A->copy(K);
        solver->A->updated = true;
        solver->b->copy(f);
        solver->b->add(-1, R);
        solver->dirichlet->copy(dirichlet);
        solver->dirichlet->add(-1, U);

        eslog::info("      == ----------------------------------------------------------------------------- == \n");
        eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s == \n", eslog::time() - start);

        solver->update(step);
        solver->solve(step);

        if ((converged = checkTemp(step))) {
            break;
        }
    }
    info::mesh->output->updateSolution(step, time);
    return converged;
}

bool HeatSteadyStateNonLinear::checkTemp(step::Step &step)
{
    double solution = eslog::time();
    double solutionNumerator = solver->x->norm();
    solver->x->add(1, U);
    x->copy(solver->x);

    double solutionDenominator = std::max(solver->x->norm(), 1e-3);
    double norm = solutionNumerator / solutionDenominator;

    eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
    eslog::info("      == ----------------------------------------------------------------------------- == \n");

    if (norm > configuration.nonlinear_solver.requested_first_residual) {
        eslog::info("      == TEMPERATURE NORM, CRITERIA                          %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.nonlinear_solver.requested_first_residual);
    } else {
        eslog::info("      == TEMPERATURE NORM, CRITERIA [CONVERGED]              %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.nonlinear_solver.requested_first_residual);
    }

    assembler.updateSolution(step, x);
    return !(norm > configuration.nonlinear_solver.requested_first_residual);
}

void HeatSteadyStateNonLinear::storeSystem(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{K, R, f}\n");
        K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
        R->store(utils::filename(utils::debugDirectory(step) + "/scheme", "R").c_str());
        f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
        dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet").c_str());
    }
}

void HeatSteadyStateNonLinear::storeSolution(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{x}\n");
        x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
    }
}


