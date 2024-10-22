
#include "structuralmechanics.transient.linear.h"
#include "structuralmechanics.steadystate.linear.h"
#include "structuralmechanics.steadystate.nonlinear.h"

#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/surfacestore.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"
#include "wrappers/precice/w.precice.h"

using namespace espreso;

StructuralMechanicsTransientLinear::StructuralMechanicsTransientLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration},
  K{}, M{}, f{}, x{}, checkpoint{},
  U{}, dU{}, V{}, W{}, X{}, Y{}, Z{}, dTK{}, dTM{},
  dirichlet{},
  pattern{}, solver{}
{

}

StructuralMechanicsTransientLinear::~StructuralMechanicsTransientLinear()
{
    if (K) { delete K; }
    if (M) { delete M; }
    if (f) { delete f; }
    if (x) { delete x; }
    if (dirichlet) { delete dirichlet; }
    if (checkpoint) { delete checkpoint; }

    if (  U) { delete   U; }
    if ( dU) { delete  dU; }
    if (  V) { delete   V; }
    if (  W) { delete   W; }
    if (  X) { delete   X; }
    if (  Y) { delete   Y; }
    if (  Z) { delete   Z; }
    if (dTK) { delete dTK; }
    if (dTM) { delete dTM; }

    if (pattern) { delete pattern; }
    if (solver) { delete solver; }
}

bool StructuralMechanicsTransientLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == ANALYSIS                                                                      TRANSIENT == \n");
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
    M = solver->A->copyPattern();
    f = solver->b->copyPattern();
    x = solver->x->copyPattern();
    checkpoint = solver->x->copyPattern();
    dirichlet = solver->dirichlet->copyPattern();

      U = solver->b->copyPattern();
     dU = solver->b->copyPattern();
      V = solver->b->copyPattern();
      W = solver->b->copyPattern();
      X = solver->b->copyPattern();
      Y = solver->b->copyPattern();
      Z = solver->b->copyPattern();
    dTK = solver->b->copyPattern();
    dTM = solver->b->copyPattern();

    pattern->map(K);
    pattern->map(M);
    pattern->map(f);
    pattern->map(dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
    return true;
}

bool StructuralMechanicsTransientLinear::run(step::Step &step, Physics *prev)
{
    Precice precice;

    time.start = time.previous = time.current = 0;
    U->set(0);
    V->set(0);
    W->set(0);
    if (prev) {
        bool correct = false;
        if (dynamic_cast<StructuralMechanicsTransientLinear*>(prev)) {
            correct = true;
            StructuralMechanicsTransientLinear* _prev = dynamic_cast<StructuralMechanicsTransientLinear*>(prev);
            time.start = time.previous = time.current = _prev->time.final;
            U->copy(_prev->U);
            V->copy(_prev->V);
            W->copy(_prev->W);
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

    assembler.connect(K, M, f, nullptr, dirichlet);

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

    dU->set(0);
    Z->set(0);
    bool solved = true;
    while (solved && time.current + time.shift <= time.final + time.precision) {
        time.shift = precice.timeStep(time.shift);
        time.current = time.previous + time.shift;
        double newmark[] = {
            1. / (delta * time.shift * time.shift),
            alpha / (delta * time.shift),
            1. / (delta * time.shift),
            1. / (2 * delta) - 1,
            alpha / delta - 1,
            time.shift / 2 * (alpha / delta - 2),
            time.shift * (1 - alpha),
            time.shift * alpha
        };

        eslog::info(" ============================================================================================= \n");
        eslog::info(" = LOAD STEP %2d, SUBSTEP   %3d,   TIME %9.4f, TIME SHIFT %9.4f, FINAL TIME %9.4f = \n", step.loadstep + 1, step.substep + 1, time.current, time.shift, time.final);
        eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
        double start = eslog::time();
        if (precice.requiresWritingCheckpoint()) {
            checkpoint->copy(x);
        }
        precice.read(StructuralMechanics::Results::fluidForce->data.data(), time.shift);
        assembler.evaluate(step, time, K, M, f, nullptr, dirichlet);
        eslog::checkpointln("SIMULATION: PHYSICS ASSEMBLED");

        if (K->updated || M->updated || solver->A->updated) {
            solver->A->set(0)->add(1, K)->add(newmark[0], M);
            solver->A->updated = true;
        }

        X->set(0);
        X->add(newmark[0], U);
        X->add(newmark[2], V);
        X->add(newmark[3], W);
        M->apply(1, X, 0, Y);

        solver->b->copy(f);
        solver->b->add(1, Y);
        solver->b->updated = true;

        solver->dirichlet->copy(dirichlet);
        solver->dirichlet->updated = dirichlet->updated;

        storeSystem(step);

        eslog::info("       = ----------------------------------------------------------------------------- = \n");
        eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

        solver->update(step);
        eslog::checkpointln("SIMULATION: LINEAR SYSTEM UPDATED");
        solved = solver->solve(step);
        eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");
        double solution = eslog::time();

        x->copy(solver->x);
        storeSolution(step);
        assembler.updateSolution(x);

        precice.write(StructuralMechanics::Results::displacement->data.data());
        precice.advance(time.shift);
        if (precice.requiresReadingCheckpoint()) {
            x->copy(checkpoint);

            eslog::info("       = TIME STEP RESTARTED                                                %8.3f s = \n", eslog::time() - solution);
            eslog::info("       = ----------------------------------------------------------------------------- = \n");
            eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
            eslog::checkpointln("SIMULATION: SOLUTION RESTARTED");
            time.current = time.previous;
        } else {
            info::mesh->output->updateSolution(step, time);

            dU->copy(solver->x);
            dU->add(-1, U);
            U->copy(solver->x);
            Z->set(0)->add(newmark[0], dU)->add(-newmark[2], V)->add(-newmark[3], W);
            V->add(newmark[6], W)->add(newmark[7], Z);
            W->copy(Z);

            eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
            eslog::info("       = ----------------------------------------------------------------------------- = \n");
            eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
            eslog::checkpointln("SIMULATION: SOLUTION PROCESSED");
            time.previous = time.current;
            ++step.substep;
        }
    }
    return solved;
}

void StructuralMechanicsTransientLinear::storeSystem(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{K, f}\n");
        K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
        M->store(utils::filename(utils::debugDirectory(step) + "/scheme", "M").c_str());
        f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
        dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet").c_str());

        X->store(utils::filename(utils::debugDirectory(step) + "/scheme", "X").c_str());
        Y->store(utils::filename(utils::debugDirectory(step) + "/scheme", "Y").c_str());
        Z->store(utils::filename(utils::debugDirectory(step) + "/scheme", "Z").c_str());
        U->store(utils::filename(utils::debugDirectory(step) + "/scheme", "U").c_str());
        V->store(utils::filename(utils::debugDirectory(step) + "/scheme", "V").c_str());
        W->store(utils::filename(utils::debugDirectory(step) + "/scheme", "W").c_str());
        dU->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dU").c_str());
    }
}

void StructuralMechanicsTransientLinear::storeSolution(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{x}\n");
        x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
    }
}
