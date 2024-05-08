
#include "physics.h"
#include "acoustic.real.linear.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "esinfo/systeminfo.h"
#include "config/ecf/physics/acoustic.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

AcousticRealLinear::AcousticRealLinear(AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, M{}, C{}, re{}, im{}
{

}

void AcousticRealLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == ANALYSIS                                                                  HARMONIC REAL == \n");
    eslog::info(" == PHYSICS                                                                        ACOUSTIC == \n");
    eslog::info(" ============================================================================================= \n");

    step.type = step::TYPE::FREQUENCY;
    assembler.analyze(step);
    info::mesh->output->updateMonitors(step);
}

void AcousticRealLinear::run(step::Step &step)
{
//    switch (configuration.solver) {
//    case LoadStepSolverConfiguration::SOLVER::FETI:    system = new FETISystem<AcousticRealLinear>(this); break;
//    case LoadStepSolverConfiguration::SOLVER::HYPRE:   break;
//    case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new MKLPDSSSystem<AcousticRealLinear>(this); break;
//    case LoadStepSolverConfiguration::SOLVER::PARDISO: break;
//    case LoadStepSolverConfiguration::SOLVER::SUPERLU: break;
//    case LoadStepSolverConfiguration::SOLVER::WSMP:    break;
//    }

    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
//    system->setMapping(K = system->assembler.A->copyPattern());
//    system->setMapping(M = system->assembler.A->copyPattern());
//    system->setMapping(C = system->assembler.A->copyPattern());
//    system->setMapping(re.f = system->assembler.b->copyPattern());
//    system->setMapping(im.f = system->assembler.b->copyPattern());
//    system->setMapping(re.x = system->assembler.b->copyPattern());
//    system->setMapping(im.x = system->assembler.b->copyPattern());
//    system->setDirichletMapping(re.dirichlet = system->assembler.dirichlet->copyPattern());
//    system->setDirichletMapping(im.dirichlet = system->assembler.dirichlet->copyPattern());

    assembler.connect(K, M, C, re.f, im.f, nullptr, nullptr, re.dirichlet);
    frequency.shift = (configuration.harmonic_solver.max_frequency - configuration.harmonic_solver.min_frequency) / configuration.harmonic_solver.num_samples;
    frequency.start = configuration.harmonic_solver.min_frequency;
    frequency.current = frequency.start;
    frequency.angular = 2 * M_PI * frequency.current;
    frequency.final = configuration.harmonic_solver.max_frequency;

    if (MPITools::node->rank == 0) {
        info::system::memory::physics = info::system::memoryAvail();
    }
    eslog::info("  PHYSICAL SOLVER MEMORY FOOTPRINT [GB] %53.2f  \n", (info::system::memory::mesh - info::system::memory::physics) / 1024. / 1024.);
    eslog::info(" ============================================================================================= \n");

    eslog::info("\n ============================================================================================= \n");
    eslog::info(" = RUN THE SOLVER                       FREQUENCY: MIN %10.4f, MAX %10.4f, STEPS %3d = \n", configuration.harmonic_solver.min_frequency, configuration.harmonic_solver.max_frequency, configuration.harmonic_solver.num_samples);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
//    system->set(step);
    eslog::info(" ============================================================================================= \n\n");

    while (frequency.current != frequency.final) {
        double start = eslog::time();
        frequency.current += frequency.shift;
        if (frequency.current + frequency.precision >= frequency.final) {
            frequency.current = frequency.final;
        }
        frequency.angular = 2 * M_PI * frequency.current;
        eslog::info(" ============================================================================================= \n");
//        eslog::info(" = LOAD STEP %2d                                                         FREQUENCY %10.4f = \n", step::step.loadstep + 1, frequency.current);
        eslog::info(" = ----------------------------------------------------------------------------------------- = \n");

        assembler.evaluate(step, frequency, K, M, C, re.f, im.f, nullptr, nullptr, re.dirichlet);
        storeSystem(step);

        // A = [
        // K - omega^2 * M                -iC
        //              iC    K - omega^2 * M
        // ]
//        system->solver.A->touched = true;
//        system->solver.A->set(0);
//        system->solver.A->copySliced(K, 0, 0, 1, 2);
//        system->solver.A->copySliced(K, 1, 1, 1, 2);
//        system->solver.A->addSliced(-frequency.angular * frequency.angular, M, 0, 0, 1, 2);
//        system->solver.A->addSliced(-frequency.angular * frequency.angular, M, 1, 1, 1, 2);
//
//        system->solver.A->addSliced(-frequency.angular, C, 0, 1, 1, 2);
//        system->solver.A->addSliced( frequency.angular, C, 1, 0, 1, 2);
//
//        system->solver.b->touched = true;
//        system->solver.b->set(0);
//        system->solver.b->copySliced(re.f, 0, 1, 2);
//        system->solver.b->copySliced(im.f, 1, 1, 2);
//
//        system->solver.dirichlet->touched = true;
//        system->solver.dirichlet->set(0);
//        system->solver.dirichlet->copySliced(re.dirichlet, 0, 1, 2);
//        system->solver.dirichlet->copySliced(im.dirichlet, 1, 1, 2);

        eslog::info("       = ----------------------------------------------------------------------------- = \n");
        eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);
//        system->update(step);
//        system->solve(step);

        double solution = eslog::time();
//        re.x->copySliced(system->solver.x, 0, 1, 2);
//        im.x->copySliced(system->solver.x, 1, 1, 2);
        assembler.updateSolution(step, re.x, im.x);
        info::mesh->output->updateSolution(step, frequency);
        eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
        eslog::info("       = ----------------------------------------------------------------------------- = \n");

        eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
    }
}

void AcousticRealLinear::storeSystem(step::Step &step)
{
//    if (info::ecf->output.print_matrices) {
//        eslog::storedata(" STORE: scheme/{K, M, C, f.re, f.im, dirichlet.re, dirichlet.im}\n");
//        K->store(utils::filename(utils::debugDirectory() + "/scheme", "K").c_str());
//        M->store(utils::filename(utils::debugDirectory() + "/scheme", "M").c_str());
//        C->store(utils::filename(utils::debugDirectory() + "/scheme", "C").c_str());
//        re.f->store(utils::filename(utils::debugDirectory() + "/scheme", "f.re").c_str());
//        im.f->store(utils::filename(utils::debugDirectory() + "/scheme", "f.im").c_str());
//        re.dirichlet->store(utils::filename(utils::debugDirectory() + "/scheme", "dirichlet.re").c_str());
//        im.dirichlet->store(utils::filename(utils::debugDirectory() + "/scheme", "dirichlet.im").c_str());
//    }
}

void AcousticRealLinear::storeSolution(step::Step &step)
{
//    if (info::ecf->output.print_matrices) {
//        eslog::storedata(" STORE: scheme/{x.re, x.im}\n");
//        re.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.re").c_str());
//        im.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.im").c_str());
//    }
}
