

#include "structuralmechanics.transient.nonlinear.h"
#include "structuralmechanics.steadystate.linear.h"
#include "structuralmechanics.steadystate.nonlinear.h"

#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"
#include "wrappers/precice/w.precice.h"

using namespace espreso;

StructuralMechanicsTransientNonLinear::StructuralMechanicsTransientNonLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration},
  K{}, M{}, C{}, f{}, f_old{},
  R{}, R_old{}, dU{}, U{}, V{}, A{}, U_old{}, V_old{}, A_old{}, X{},
  dirichlet{},
  pattern{}, solver{}
{

}

StructuralMechanicsTransientNonLinear::~StructuralMechanicsTransientNonLinear()
{
    if (K) { delete K; }
    if (M) { delete M; }
    if (C) { delete C; }
    if (f) { delete f; }
    if (f_old) { delete f_old; }
    if (dirichlet) { delete dirichlet; }

    if (  R) { delete   R; }
    if (R_old) { delete R_old; }
    if ( dU) { delete  dU; }
    if (  U) { delete   U; }
    if (  V) { delete   V; }
    if (  A) { delete   A; }
    if (U_old) { delete U_old; }
    if (V_old) { delete V_old; }
    if (A_old) { delete A_old; }
    if (X) { delete X; }

    if (pattern) { delete pattern; }
    if (solver) { delete solver; }
}

bool StructuralMechanicsTransientNonLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == ANALYSIS                                                                      TRANSIENT == \n");
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
    M = solver->A->copyPattern();
    C = solver->A->copyPattern();
    f = solver->b->copyPattern();
    f_old = solver->b->copyPattern();
    dirichlet = solver->dirichlet->copyPattern();

    R = solver->b->copyPattern();
    R_old = solver->b->copyPattern();
    dU = solver->b->copyPattern();
    U = solver->b->copyPattern();
    V = solver->b->copyPattern();
    A = solver->b->copyPattern();
    U_old = solver->b->copyPattern();
    V_old = solver->b->copyPattern();
    A_old = solver->b->copyPattern();
    X = solver->b->copyPattern();

    pattern->map(K);
    pattern->map(M);
    pattern->map(C);
    pattern->map(f);
    pattern->map(R);
    pattern->map(dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
    return true;
}

bool StructuralMechanicsTransientNonLinear::run(step::Step &step, Physics *prev)
{
    Precice precice;

    time.start = time.previous = time.current = 0;
    U->set(0);
    V->set(0);
    A->set(0);
    if (prev) {
        bool correct = false;
        if (dynamic_cast<StructuralMechanicsTransientNonLinear*>(prev)) {
            correct = true;
            StructuralMechanicsTransientNonLinear* _prev = dynamic_cast<StructuralMechanicsTransientNonLinear*>(prev);
            time.start = time.previous = time.current = _prev->time.final;
            U->copy(_prev->U);
            V->copy(_prev->V);
            A->copy(_prev->A);
        }
        if (dynamic_cast<StructuralMechanicsSteadyStateLinear*>(prev)) {
            correct = true;
            StructuralMechanicsSteadyStateLinear* _prev = dynamic_cast<StructuralMechanicsSteadyStateLinear*>(prev);
            time.start = time.previous = time.current = _prev->time.final;
            U->copy(_prev->x);
        }
        if (dynamic_cast<StructuralMechanicsSteadyStateNonLinear*>(prev)) {
            correct = true;
            StructuralMechanicsSteadyStateNonLinear* _prev = dynamic_cast<StructuralMechanicsSteadyStateNonLinear*>(prev);
            time.start = time.previous = time.current = _prev->time.final;
            U->copy(_prev->x);
        }
        if (!correct) {
            eslog::globalerror("Incompatible load steps.\n");
        }
        assembler.updateSolution(U);
    } else {
        assembler.getInitialVelocity(V);
    }

    time.shift = configuration.transient_solver.time_step;
    time.final = time.start + configuration.duration_time;

    double alpha = configuration.transient_solver.alpha;
    double delta = configuration.transient_solver.delta;
    double alphaM = configuration.transient_solver.alphaM;
    double alphaF = configuration.transient_solver.alphaF;

    double dampingK = configuration.transient_solver.damping.rayleigh.direct_damping.stiffness.evaluator->evaluate();
    double dampingM = configuration.transient_solver.damping.rayleigh.direct_damping.mass.evaluator->evaluate();

    double nd = configuration.transient_solver.numerical_damping;
    alpha *= (1 + nd) * (1 + nd);
    delta += nd;

    assembler.connect(K, M, f, R, dirichlet);

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
    eslog::info(" = TRANSIENT SOLVER CONFIGURATION                                                            = \n");
    eslog::info(" =  - ALPHA                                                                        %9.4f = \n", alpha);
    eslog::info(" =  - DELTA                                                                        %9.4f = \n", delta);
    eslog::info(" =  - ALPHA M                                                                      %9.4f = \n", alphaM);
    eslog::info(" =  - ALPHA F                                                                      %9.4f = \n", alphaF);
    eslog::info(" =  - DIRECT STIFFNESS DAMPING                                                     %9.4f = \n", dampingK);
    eslog::info(" =  - DIRECT MASS DAMPING                                                          %9.4f = \n", dampingM);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    eslog::info(" = NEWTON RAPHSON CONVERGENCE CRITERIA                                                       = \n");
    if (configuration.nonlinear_solver.check_first_residual) {
        eslog::info(" =  - DISPLACEMENT RESIDUAL                                                             TRUE = \n");
        eslog::info(" =  - DISPLACEMENT RESIDUAL PRECISION                                           %e = \n", configuration.nonlinear_solver.requested_first_residual);
    } else {
        eslog::info(" =  - DISPLACEMENT RESIDUAL                                                            FALSE = \n");
    }
    if (configuration.nonlinear_solver.check_second_residual) {
        eslog::info(" =  - STRESS RESIDUAL                                                                   TRUE = \n");
        eslog::info(" =  - STRESS RESIDUAL PRECISION                                                 %e = \n", configuration.nonlinear_solver.requested_second_residual);
    } else {
        eslog::info(" =  - STRESS RESIDUAL                                                                  FALSE = \n");
    }

    eslog::info(" ============================================================================================= \n");
    eslog::info(" = LOAD STEP %2d, INITIALIZATION   TIME %9.4f, TIME SHIFT %9.4f, FINAL TIME %9.4f = \n", step.loadstep + 1, time.current, time.shift, time.final);
    eslog::info(" ============================================================================================= \n");

    double tinit = eslog::time();
    C->set(0);
    dU->set(0);
    solver->set(step);

    assembler.evaluate(step, time, K, M, f, nullptr, dirichlet);
    C->add(dampingK, K);
    C->add(dampingM, M);
    f_old->copy(f);
    U_old->copy(U);
    V_old->copy(V);
    A_old->copy(A);

    eslog::info(" = INITIALIZATION                                                                 %8.3f s = \n", eslog::time() - tinit);
    eslog::checkpointln("SIMULATION: TRANSIENT SOLVER INITIALIZED");

    step.substep = 0;
    bool converged = true;
    while (converged && time.current + time.shift <= time.final + time.precision) {
        step.iteration = 0;
        double start = eslog::time();
        time.shift = precice.timeStep(time.shift);
        time.current = time.previous + time.shift;

        eslog::info(" ============================================================================================= \n");
        eslog::info(" == SUBSTEP   %3d               TIME %8.5f    TIME SHIFT %8.5f    FINAL TIME %8.5f == \n", step.substep, time.current, time.shift, time.final);
        eslog::info(" ============================================================================================= \n");

        double a0 = (1. - alphaM) / (alpha * time.shift * time.shift);
        double a1 = ((1 - alphaF) * time.shift) / (alpha * time.shift);
        double a2 = a0 * time.shift;
        double a3 = (1 - alphaM) / (2 * alpha) - 1;
        double a4 = ((1 - alphaF) * time.shift) / alpha - 1;
        double a5 = (1 - alphaF) * (time.shift / (2 * alpha) - 1) * time.shift;

        if (precice.requiresWritingCheckpoint()) {
            // previous data are stored in *_old
        }
        precice.read(StructuralMechanics::Results::fluidForce->data.data(), time.shift);

        assembler.evaluate(step, time, K, nullptr, f, R, dirichlet);
        eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

        solver->A->set(0)->add(a0, M)->add(a1, C)->add(1 - alphaF, K);
        solver->A->updated = true;

        solver->b->set(0)->add(1 - alphaF, f)->add(alphaF, f_old)->add(-1., R);
        X->set(0)->add(a2, V_old)->add(a3, A_old); M->apply(1., X, 1., solver->b);
        X->set(0)->add(a4, V_old)->add(a5, A_old); C->apply(1., X, 1., solver->b);
        solver->b->updated = true;

        solver->dirichlet->copy(dirichlet);
        solver->dirichlet->add(-1, U_old);
        solver->dirichlet->updated = true;

        storeSystem(step);
        solver->update(step);
        eslog::checkpointln("SIMULATION: LINEAR SYSTEM UPDATED");
        solver->solve(step);
        storeSolution(step);
        eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");

        dU->copy(solver->x);
        U->copy(U_old)->add(1., dU);
        V->set(0);
        V->add(delta / (alpha * time.shift), dU);
        V->add(-delta / alpha + 1, V_old);
        V->add(-time.shift / 2 * (delta / alpha - 2), A_old);
        A->set(0);
        A->add(1. / (alpha * time.shift * time.shift), dU);
        A->add(-1. / (alpha * time.shift), V_old);
        A->add(-1. / (2 * alpha) - 1, A_old);
        R_old->copy(R);

        assembler.updateSolution(U);
        eslog::checkpointln("SIMULATION: SOLUTION UPDATED");
        eslog::info("       = PREDICTOR COMPUTED                                                 %8.3f s = \n", eslog::time() - start);
        eslog::info("       = ----------------------------------------------------------------------------- = \n");

        // NEWTON RAPHSON
        converged = false;
        double f_norm = f->norm();
        double U_norm = dU->norm(); // U from the predictor step
        while (!converged && step.iteration++ < configuration.nonlinear_solver.max_iterations) {
            double istart = eslog::time();
            eslog::info("      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

            assembler.evaluate(step, time, K, nullptr, nullptr, R, nullptr);
            eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - istart);

            // KEffCor = (a0 * M + a1 * C + (1 - alphaF) * Kt_new);
            solver->A->set(0)->add(a0, M)->add(a1, C)->add(1 - alphaF, K);
            solver->A->updated = true;

            // rEffCor = (1 - alphaF) * (f_ext_new - f_int_new - C * v_new) + alphaF * (f_ext_old - f_int_old - C * v_old) - M * ((1 - alphaM) * a_new + alphaM * a_old);
            solver->b->set(0);
            solver->b->add(1 - alphaF, f)->add(alphaF - 1, R); C->apply(alphaF - 1, V, 1., solver->b);
            solver->b->add(alphaF, f_old)->add(-alphaF, R_old); C->apply(-alphaF, V_old, 1., solver->b);
            X->set(0)->add(1 - alphaM, A)->add(alphaM, A_old); M->apply(-1., X, 1., solver->b);
            solver->b->updated = true;

            storeSystem(step);
            solver->update(step);
            eslog::checkpointln("SIMULATION: LINEAR SYSTEM UPDATED");
            solver->solve(step);
            storeSolution(step);
            eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");

            dU->copy(solver->x);
            U->add(1., dU);
            V->set(0);
            V->add(delta / (alpha * time.shift), U)->add(-delta / (alpha * time.shift), U_old);
            V->add(-delta / alpha + 1, V_old)->add(-time.shift / 2 * (delta / alpha - 2), A_old);
            A->set(0);
            A->add(1. / (alpha * time.shift * time.shift), U)->add(-1. / (alpha * time.shift * time.shift), U_old);
            A->add(-1. / (alpha * time.shift), V_old)->add(-1. / (2 * alpha) + 1, A_old);
            eslog::checkpointln("SIMULATION: SOLUTION UPDATED");

            bool norm1 = checkDisplacement(step, U_norm);
            bool norm2 = checkStress(step, f_norm);
            converged = norm1 && norm2;
            if (converged) {
                eslog::info("      =================================================================================== \n\n");
                assembler.updateSolution(U);
            } else {
                eslog::info("       = ----------------------------------------------------------------------------- = \n");
                assembler.nextIteration(U);
            }
        }
        precice.write(StructuralMechanics::Results::displacement->data.data());
        precice.advance(time.shift);
        if (precice.requiresReadingCheckpoint()) {
            assembler.updateSolution(U_old);
            eslog::info("      = TIME STEP RESTARTED ====================================== solved in %8.3f s = \n", eslog::time() - start);
            eslog::info("      =================================================================================== \n\n");
            eslog::checkpointln("SIMULATION: SOLUTION RESTARTED");
            time.current = time.previous;
        } else {
            info::mesh->output->updateSolution(step, time);

            U_old->copy(U);
            V_old->copy(V);
            A_old->copy(A);
            f_old->copy(f);

            if (StructuralMechanics::Results::acceleration) {
                A->storeTo(StructuralMechanics::Results::acceleration->data);
            }
            if (StructuralMechanics::Results::velocity) {
                V->storeTo(StructuralMechanics::Results::velocity->data);
            }

            eslog::info("      = TIME STEP FINISHED ======================================= solved in %8.3f s = \n", eslog::time() - start);
            eslog::info("      =================================================================================== \n\n");
            eslog::checkpointln("SIMULATION: SOLUTION PROCESSED");
            time.previous = time.current;
            ++step.substep;
        }
    }
    return converged;
}

bool StructuralMechanicsTransientNonLinear::checkDisplacement(step::Step &step, double U_norm)
{
    if (!configuration.nonlinear_solver.check_first_residual) {
        return true;
    }
    double solutionNumerator = solver->x->norm();
    double solutionDenominator = std::max(U_norm, 1e-3);
    double norm = solutionNumerator / solutionDenominator;

    if (norm > configuration.nonlinear_solver.requested_first_residual) {
        eslog::info("      == DISPLACEMENT NORM / CRITERIA                        %.5e / %.5e == \n", norm, configuration.nonlinear_solver.requested_first_residual);
        return false;
    } else {
        eslog::info("      == DISPLACEMENT NORM / CRITERIA [CONVERGED]            %.5e / %.5e == \n", norm, configuration.nonlinear_solver.requested_first_residual);
        return true;
    }
}

bool StructuralMechanicsTransientNonLinear::checkStress(step::Step &step, double f_norm)
{
    if (!configuration.nonlinear_solver.check_second_residual) {
        return true;
    }
    double b_norm = solver->rhs_norm();
    double nR = b_norm / (1 + f_norm);

    if (nR > configuration.nonlinear_solver.requested_second_residual) {
        eslog::info("      == STRESS NORM       / CRITERIA                        %.5e / %.5e == \n", nR, configuration.nonlinear_solver.requested_second_residual);
        return false;
    } else {
        eslog::info("      == STRESS NORM       / CRITERIA [CONVERGED]            %.5e / %.5e == \n", nR, configuration.nonlinear_solver.requested_second_residual);
        return true;
    }
}

void StructuralMechanicsTransientNonLinear::storeSystem(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{K, f}\n");
        K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
        M->store(utils::filename(utils::debugDirectory(step) + "/scheme", "M").c_str());
        C->store(utils::filename(utils::debugDirectory(step) + "/scheme", "C").c_str());
        f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
        dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet").c_str());

        R->store(utils::filename(utils::debugDirectory(step) + "/scheme", "R").c_str());
        dU->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dU").c_str());
        U->store(utils::filename(utils::debugDirectory(step) + "/scheme", "U").c_str());
        V->store(utils::filename(utils::debugDirectory(step) + "/scheme", "V").c_str());
        A->store(utils::filename(utils::debugDirectory(step) + "/scheme", "A").c_str());
    }
}

void StructuralMechanicsTransientNonLinear::storeSolution(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{x}\n");
        solver->x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
    }
}
