

#include "structuralmechanics.transient.nonlinear.h"

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
#include "wrappers/precice/w.precice.h"

using namespace espreso;

StructuralMechanicsTransientNonLinear::StructuralMechanicsTransientNonLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration},
  K{}, M{}, C{}, f{}, f_old{}, x{}, dirichlet{}, prev{},
  R{}, R_old{}, dU{}, U{}, V{}, A{}, U_old{}, V_old{}, A_old{}, X{},
  builder{}, solver{}
{

}

StructuralMechanicsTransientNonLinear::~StructuralMechanicsTransientNonLinear()
{
    if (K) { delete K; }
    if (M) { delete M; }
    if (C) { delete C; }
    if (f) { delete f; }
    if (f_old) { delete f_old; }
    if (x) { delete x; }
    if (dirichlet) { delete dirichlet; }
    if (prev) { delete prev; }

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

    if (builder) { delete builder; }
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

    switch (configuration.solver) {
    case LoadStepSolverConfiguration::SOLVER::FETI:
        builder = new UniformBuilderFETI<double>(configuration, 1);
        solver = new FETILinearSystemSolver<double>(settings, configuration);
        break;
    case LoadStepSolverConfiguration::SOLVER::HYPRE:   break;
    case LoadStepSolverConfiguration::SOLVER::MKLPDSS:
        builder = new UniformBuilderDirect<double>(configuration, 1);
        solver = new MKLPDSSLinearSystemSolver<double>(configuration.mklpdss);
        break;
    case LoadStepSolverConfiguration::SOLVER::PARDISO: break;
    case LoadStepSolverConfiguration::SOLVER::SUPERLU: break;
    case LoadStepSolverConfiguration::SOLVER::WSMP:    break;
    case LoadStepSolverConfiguration::SOLVER::NONE:
        builder = new UniformBuilderDirect<double>(configuration, 1);
        solver = new EmptySystemSolver<double>();
    }

    builder->fillMatrix(solver->A);
    builder->fillVector(solver->b);
    builder->fillVector(solver->x);
    builder->fillDirichlet(solver->dirichlet);

    K = solver->A->copyPattern();
    M = solver->A->copyPattern();
    C = solver->A->copyPattern();
    f = solver->b->copyPattern();
    x = solver->x->copyPattern();
    prev = solver->x->copyPattern();
    dirichlet = solver->dirichlet->copyPattern();

      R = solver->b->copyPattern();
     dU = solver->b->copyPattern();
      U = solver->b->copyPattern();
      V = solver->b->copyPattern();
      A = solver->b->copyPattern();

    builder->fillMatrixMap(K);
    builder->fillMatrixMap(M);
    builder->fillMatrixMap(C);
    builder->fillVectorMap(f);
    builder->fillVectorMap(R);
    builder->fillDirichletMap(dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
    return true;
}

bool StructuralMechanicsTransientNonLinear::run(step::Step &step)
{
    Precice precice;

    time.start = time.previous = time.current = 0;
    time.shift = configuration.transient_solver.time_step;
    time.final = configuration.duration_time;

    double alpha = configuration.transient_solver.alpha;
    double delta = configuration.transient_solver.delta;
    double alphaM = configuration.transient_solver.alphaM;
    double alphaF = configuration.transient_solver.alphaF;

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
    eslog::info(" = NEWTON RAPHSON CONVERGENCE CRITERIA                                                      = \n");
    if (configuration.nonlinear_solver.check_first_residual) {
        eslog::info(" =  - DISPLACEMENT RESIDUAL                                                             TRUE = \n");
        eslog::info(" =  - DISPLACEMENT RESIDUAL PRECISION                                                     %e = \n", configuration.nonlinear_solver.requested_first_residual);
    } else {
        eslog::info(" =  - DISPLACEMENT RESIDUAL                                                            FALSE = \n");
    }
    if (configuration.nonlinear_solver.check_second_residual) {
        eslog::info(" =  - STRESS RESIDUAL                                                                   TRUE = \n");
        eslog::info(" =  - STRESS RESIDUAL PRECISION                                                           %e = \n", configuration.nonlinear_solver.requested_second_residual);
    } else {
        eslog::info(" =  - STRESS RESIDUAL                                                                  FALSE = \n");
    }

    eslog::info(" ============================================================================================= \n");
    eslog::info(" = LOAD STEP %2d, INIT              TIME %9.4f, TIME SHIFT %9.4f, FINAL TIME %9.4f = \n", step.loadstep + 1, time.current, time.shift, time.final);
    eslog::info("      =================================================================================== \n\n");

    C->set(0);
    dU->set(0);
    U->set(0);
    V->set(0);
    A->set(0);
    solver->set(step);

    assembler.evaluate(step, time, K, M, f, nullptr, dirichlet);
    C->add(configuration.transient_solver.damping.rayleigh.direct_damping.stiffness.evaluator->evaluate(), K);
    C->add(configuration.transient_solver.damping.rayleigh.direct_damping.mass.evaluator->evaluate(), M);
    f_old->copy(f);
    U_old->copy(U);
    V_old->copy(V);
    A_old->copy(A);

    step.substep = 0;
    bool converged = true;
    while (converged && time.current + time.shift <= time.final + time.precision) {
        double start = eslog::time();
        time.shift = precice.timeStep(time.shift);
        time.current = time.previous + time.shift;

        double a0 = (1. - alphaM) / (alpha * time.shift * time.shift);
        double a1 = ((1 - alphaF) * time.shift) / (alpha * time.shift);
        double a2 = a0 * time.shift;
        double a3 = (1 - alphaM) / (2 * alpha) - 1;
        double a4 = ((1 - alphaF) * time.shift) / alpha - 1;
        double a5 = (1 - alphaF) * (time.shift / (2 * alpha) - 1) * time.shift;


        if (precice.requiresWritingCheckpoint()) {
//            prev->copy(U);
        }
        precice.read(StructuralMechanics::Results::fluidForce->data.data(), time.shift);

        eslog::info("\n      ==                                                        INITIAL ITERATION == \n", step.iteration);
        assembler.evaluate(step, time, K, nullptr, f, R, dirichlet);

        solver->A->set(0)->add(a0, M)->add(a1, C)->add(1 - alphaF, K);
        solver->A->updated = true;

        solver->b->set(0)->add(1 - alphaF, f)->add(alphaF, f_old)->add(-1., R);
        X->set(0)->add(a2, V_old)->add(a3, A_old); M->apply(1., X, 1., solver->b);
        X->set(0)->add(a4, V_old)->add(a5, A_old); C->apply(1., X, 1., solver->b);
        solver->b->updated = true;

        solver->dirichlet->copy(dirichlet);
        solver->dirichlet->add(-1, U);
        solver->dirichlet->updated = true;

        storeSystem(step);
        solver->update(step);
        solver->solve(step);

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

        // NEWTON RAPHSON
        step.iteration = 0;
        converged = false;
        double f_norm = f->norm();
        while (step.iteration++ < configuration.nonlinear_solver.max_iterations) {
            eslog::info("\n      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

            assembler.evaluate(step, time, K, nullptr, nullptr, R, nullptr);
            solver->A->set(0)->add(a0, M)->add(a1, C)->add(1 - alphaF, K);
            solver->A->updated = true;

            solver->b->set(0);
            solver->b->add(1 - alphaF, f)->add(alphaF - 1, R); C->apply(-1., V, 1., solver->b);
            solver->b->add(alphaF, f_old)->add(-alphaF, R_old); C->apply(-1., V_old, 1., solver->b);
            X->copy(A)->add(1., A_old); M->apply(-1., X, 1., solver->b);
            solver->b->updated = true;

            storeSystem(step);
            solver->update(step);
            solver->solve(step);

            dU->copy(solver->x);
            U->add(1., dU);
            V->set(0);
            V->add(delta / (alpha * time.shift), U)->add(-delta / (alpha * time.shift), U_old);
            V->add(-delta / alpha + 1, V_old)->add(-time.shift / 2 * (delta / alpha - 2), A_old);
            A->set(0);
            A->add(1. / (alpha * time.shift * time.shift), U)->add(-1. / (alpha * time.shift * time.shift), U_old);
            A->add(-1. / (alpha * time.shift), V_old)->add(-1. / (2 * alpha) - 1, A_old);

            if ((converged = checkDisplacement(step, f_norm))) {
                break;
            }
        }
        precice.write(StructuralMechanics::Results::displacement->data.data());
        precice.advance(time.shift);
        if (precice.requiresReadingCheckpoint()) {
            x->copy(prev);

            eslog::info("       = TIME STEP RESTARTED                                                %8.3f s = \n", eslog::time() - start);
            eslog::info("       = ----------------------------------------------------------------------------- = \n");
            eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
            eslog::checkpointln("SIMULATION: SOLUTION RESTARTED");
            time.current = time.previous;
        } else {
            info::mesh->output->updateSolution(step, time);
            storeSolution(step);

            U_old->copy(U);
            V_old->copy(V);
            A_old->copy(A);

            eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - start);
            eslog::info("       = ----------------------------------------------------------------------------- = \n");
            eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
            eslog::checkpointln("SIMULATION: SOLUTION PROCESSED");
            time.previous = time.current;
            ++step.substep;
        }
    }
    return converged;
}

bool StructuralMechanicsTransientNonLinear::checkDisplacement(step::Step &step, double f_norm)
{
    double b_norm = solver->b->norm();
    double nR = b_norm / (1 + f_norm);

    if (nR > configuration.nonlinear_solver.requested_first_residual) {
        eslog::info("      == DISPLACEMENT NORM, CRITERIA                         %.5e / %.5e == \n", b_norm, f_norm * configuration.nonlinear_solver.requested_first_residual);
        assembler.nextIteration(U);
        return false;
    } else {
        eslog::info("      == DISPLACEMENT NORM, CRITERIA [CONVERGED]             %.5e / %.5e == \n", b_norm, f_norm * configuration.nonlinear_solver.requested_first_residual);
        eslog::info("      =================================================================================== \n\n");
        assembler.updateSolution(U);
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
        x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
    }
}
