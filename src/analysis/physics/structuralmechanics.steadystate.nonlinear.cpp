
#include "structuralmechanics.steadystate.nonlinear.h"

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
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, U{}, R{}, f{}, dirichlet{}, pattern{}, solver{}
{

}

StructuralMechanicsSteadyStateNonLinear::~StructuralMechanicsSteadyStateNonLinear()
{
    if (K) { delete K; }
    if (U) { delete U; }
    if (R) { delete R; }
    if (f) { delete f; }
    if (dirichlet) { delete dirichlet; }
    if (pattern) { delete pattern; }
    if (solver) { delete solver; }
}

bool StructuralMechanicsSteadyStateNonLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == ANALYSIS                                                                   STEADY STATE == \n");
    eslog::info(" == PHYSICS                                                            STRUCTURAL MECHANICS == \n");
    eslog::info(" == MODE                                                                         NON-LINEAR == \n");
    eslog::info(" ============================================================================================= \n");

    step.type = step::TYPE::TIME;
    if (!assembler.analyze(step)) {
        return false;
    }
    info::mesh->output->updateMonitors(step);
    if (configuration.nonlinear_solver.stepping == NonLinearSolverConfiguration::STEPPINGG::FALSE) {
        configuration.nonlinear_solver.substeps = 1;
    }

    solver = setSolver<double>(settings, configuration);
    pattern = solver->getPattern(configuration, 1);

    pattern->set(solver->A);
    pattern->set(solver->b);
    pattern->set(solver->x);
    pattern->set(solver->dirichlet);

    K = solver->A->copyPattern();
    R = solver->b->copyPattern();
    U = solver->b->copyPattern();
    f = solver->b->copyPattern();
    dirichlet = solver->dirichlet->copyPattern();

    pattern->map(K);
    pattern->map(f);
    pattern->map(R);
    pattern->map(dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
    return true;
}

bool StructuralMechanicsSteadyStateNonLinear::run(step::Step &step, Physics *prev)
{
    if (prev) {
        bool correct = false;
        if (!correct) {
            eslog::globalerror("Incompatible load steps.\n");
        }
    }

    assembler.connect(K, nullptr, f, R, dirichlet);

    if (MPITools::node->rank == 0) {
        info::system::memory::physics = info::system::memoryAvail();
    }
    eslog::info("  PHYSICAL SOLVER MEMORY FOOTPRINT [GB] %53.2f  \n", (info::system::memory::mesh - info::system::memory::physics) / 1024. / 1024.);
    eslog::info(" ============================================================================================= \n");

    eslog::info("\n ============================================================================================= \n");
    eslog::info(" = RUN THE SOLVER                                                DURATION TIME: %10.4f s = \n", configuration.duration_time);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    eslog::info(" = MODE                                                                           NON-LINEAR = \n");
    eslog::info(" = NUMBER OF SUBSTEPS                                                             %10d = \n", configuration.nonlinear_solver.substeps);
    eslog::info(" = MAX ITERATIONS                                                                 %10d = \n", configuration.nonlinear_solver.max_iterations);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    solver->set(step);
    eslog::info(" ============================================================================================= \n\n");
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM SET");

    eslog::info(" ============================================================================================= \n");
    eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step.loadstep + 1, time.current);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    eslog::info("      =================================================================================== \n");
    eslog::info("      ==  NEWTON RAPHSON CONVERGENCE CRITERIA                                          == \n");
    if (configuration.nonlinear_solver.check_first_residual) {
        eslog::info("      ==  - DISPLACEMENT RESIDUAL                                                 TRUE == \n");
        eslog::info("      ==  - REQUESTED NORM                                                 %.5e == \n", configuration.nonlinear_solver.requested_first_residual);
    } else {
        eslog::info("      ==  - DISPLACEMENT RESIDUAL                                                FALSE == \n");
    }
    if (configuration.nonlinear_solver.check_second_residual) {
        eslog::info("      ==  - STRESS RESIDUAL                                                       TRUE == \n");
    } else {
        eslog::info("      ==  - STRESS RESIDUAL                                                      FALSE == \n");
    }
    eslog::info("      =================================================================================== \n\n");

    time.start = 0;
    time.shift = configuration.duration_time / configuration.nonlinear_solver.substeps;
    time.final = configuration.duration_time;
    step.substeps = configuration.nonlinear_solver.substeps;
    U->set(0);
    bool converged = true;
    for (step.substep = 0; converged && step.substep < configuration.nonlinear_solver.substeps; ++step.substep) {
        if (configuration.nonlinear_solver.substeps > 1) {
            eslog::info("      ==  SUBSTEP                                                           %10d ==     \n", step.substep + 1);
            eslog::info("      =================================================================================== \n");
        }
        eslog::info("      ==                                                                  INITIAL STEP ==     \n");

        double start = eslog::time();
        step.iteration = 0;
        time.current = (step.substep + 1) * time.shift;
        if (step.substep + 1 == configuration.nonlinear_solver.substeps) {
            time.current = time.final;
        }

        assembler.evaluate(step, time, K, nullptr, f, R, dirichlet);
        storeSystem(step);
        solver->A->copy(K);
        solver->A->updated = true;
        solver->b->copy(f);
        solver->dirichlet->copy(dirichlet);
        solver->dirichlet->add(-1, U);
        solver->dirichlet->updated = true;
        eslog::info("      == ----------------------------------------------------------------------------- == \n");
        eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

        solver->update(step);
        solver->solve(step);

        double solution = eslog::time();
        U->copy(solver->x);
        storeSolution(step);
        assembler.nextIteration(U);
        eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
        eslog::info("      == ----------------------------------------------------------------------------- == \n");

        converged = false;
        double f_norm = f->norm();
        double U_norm = U->norm(); // U from the linear step
        while (!converged && step.iteration++ < configuration.nonlinear_solver.max_iterations) {
            eslog::info("\n      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

            start = eslog::time();
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
            U->add(1., solver->x);

            bool norm1 = checkDisplacement(step, U_norm);
            bool norm2 = checkStress(step, f_norm);
            converged = norm1 && norm2;
            if (converged) {
                eslog::info("      =================================================================================== \n\n");
                storeSolution(step);
                assembler.updateSolution(U);
            } else {
                eslog::info("       = ----------------------------------------------------------------------------- = \n");
                storeSolution(step);
                assembler.nextIteration(U);
            }
        }
        info::mesh->output->updateSolution(step, time);
    }
    return converged;
}

bool StructuralMechanicsSteadyStateNonLinear::checkDisplacement(step::Step &step, double U_norm)
{
    if (!configuration.nonlinear_solver.check_first_residual) {
        return true;
    }
    double solutionNumerator = solver->x->norm();
    double solutionDenominator = std::max(U_norm, 1e-3);
    printf("%+.10e / %+.10e\n", solutionNumerator, U_norm);
    double norm = solutionNumerator / solutionDenominator;

    if (norm > configuration.nonlinear_solver.requested_first_residual) {
        eslog::info("      == DISPLACEMENT NORM / CRITERIA                        %.5e / %.5e == \n", norm, configuration.nonlinear_solver.requested_first_residual);
        return false;
    } else {
        eslog::info("      == DISPLACEMENT NORM / CRITERIA [CONVERGED]            %.5e / %.5e == \n", norm, configuration.nonlinear_solver.requested_first_residual);
        return true;
    }
}

bool StructuralMechanicsSteadyStateNonLinear::checkStress(step::Step &step, double f_norm)
{
    if (!configuration.nonlinear_solver.check_second_residual) {
        return true;
    }
    double b_norm = solver->rhs_without_dirichlet_norm();
    double nR = b_norm / (1 + f_norm);

    if (nR > configuration.nonlinear_solver.requested_second_residual) {
        eslog::info("      == STRESS NORM       / CRITERIA                        %.5e / %.5e == \n", nR, configuration.nonlinear_solver.requested_second_residual);
        return false;
    } else {
        eslog::info("      == STRESS NORM       / CRITERIA [CONVERGED]            %.5e / %.5e == \n", nR, configuration.nonlinear_solver.requested_second_residual);
        return true;
    }
}

void StructuralMechanicsSteadyStateNonLinear::storeSystem(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{K, R, f}\n");
        K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
        R->store(utils::filename(utils::debugDirectory(step) + "/scheme", "R").c_str());
        f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
        dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "BC").c_str());
    }
}

void StructuralMechanicsSteadyStateNonLinear::storeSolution(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{x}\n");
        U->store(utils::filename(utils::debugDirectory(step) + "/scheme", "U").c_str());
    }
}
