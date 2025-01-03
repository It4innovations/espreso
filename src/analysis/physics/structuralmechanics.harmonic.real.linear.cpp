
#include "structuralmechanics.harmonic.real.linear.h"

#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

StructuralMechanicsHarmonicRealLinear::StructuralMechanicsHarmonicRealLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration), assembler{settings, configuration}, K{}, M{}, C{}, patternAssembler{}, patternSolver{}, solver{}
{
    re.f = re.x = nullptr;
    im.f = im.x = nullptr;
    re.dirichlet = im.dirichlet = nullptr;
}

StructuralMechanicsHarmonicRealLinear::~StructuralMechanicsHarmonicRealLinear()
{
    if (K) { delete K; }
    if (M) { delete M; }
    if (C) { delete C; }
    if (re.f) { delete re.f; }
    if (im.f) { delete im.f; }
    if (re.x) { delete re.x; }
    if (im.x) { delete im.x; }
    if (re.dirichlet) { delete re.dirichlet; }
    if (im.dirichlet) { delete im.dirichlet; }
    if (patternAssembler) { delete patternAssembler; }
    if (patternSolver) { delete patternSolver; }
    if (solver) { delete solver; }
}

bool StructuralMechanicsHarmonicRealLinear::analyze(step::Step &step)
{
    eslog::info("\n ============================================================================================= \n");
    eslog::info(" == ANALYSIS                                                              HARMONIC ANALYSIS == \n");
    eslog::info(" == PHYSICS                                                            STRUCTURAL MECHANICS == \n");
    eslog::info(" ============================================================================================= \n");

    step.type = step::TYPE::FREQUENCY;
    if (!assembler.analyze(step)) {
        return false;
    }
    info::mesh->output->updateMonitors(step);

    solver = setSolver<double>(configuration);
    patternAssembler = solver->getPattern(configuration, 1);
    patternSolver = solver->getPattern(configuration, 2);

    K = solver->A->copyPattern();
    re.f = solver->b->copyPattern();
    re.dirichlet = solver->dirichlet->copyPattern();
    patternAssembler->set(K);
    patternAssembler->set(re.f);
    patternAssembler->set(re.dirichlet);
    M = K->copyPattern();
    C = K->copyPattern();
    im.f = re.f->copyPattern();
    re.x = re.f->copyPattern();
    im.x = re.f->copyPattern();
    im.dirichlet = re.dirichlet->copyPattern();

    patternSolver->set(solver);

    patternAssembler->map(K);
    patternAssembler->map(M);
    patternAssembler->map(C);
    patternAssembler->map(re.f);
    patternAssembler->map(im.f);
    patternAssembler->map(re.dirichlet);
    patternAssembler->map(im.dirichlet);
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
    return true;
}

bool StructuralMechanicsHarmonicRealLinear::run(step::Step &step, Physics *prev)
{
    step.substep = 0;
    step.substeps = configuration.harmonic_solver.num_samples;
    frequency.shift = (configuration.harmonic_solver.max_frequency - configuration.harmonic_solver.min_frequency) / configuration.harmonic_solver.num_samples;
    frequency.start = configuration.harmonic_solver.min_frequency;
    frequency.current = frequency.start;
    frequency.final = configuration.harmonic_solver.max_frequency;

    assembler.connect(K, M, C, re.f, im.f, nullptr, nullptr, re.dirichlet, im.dirichlet);

    if (MPITools::node->rank == 0) {
        info::system::memory::physics = info::system::memoryAvail();
    }
    eslog::info("  PHYSICAL SOLVER MEMORY FOOTPRINT [GB] %53.2f  \n", (info::system::memory::mesh - info::system::memory::physics) / 1024. / 1024.);
    eslog::info(" ============================================================================================= \n");

    eslog::info("\n ============================================================================================= \n");
    eslog::info(" = RUN THE SOLVER                            MIN FREQ. = %.0f, MAX FREQ. = %.0f, SAMPLES = %3d = \n", configuration.harmonic_solver.min_frequency, configuration.harmonic_solver.max_frequency, configuration.harmonic_solver.num_samples);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    solver->set(step);
    eslog::info(" ============================================================================================= \n\n");
    eslog::checkpointln("SIMULATION: LINEAR SYSTEM SET");

    int dim = info::mesh->dimension;
    Selection s1(  0, dim, 2 * dim);
    Selection s2(dim, dim, 2 * dim);
    bool solved = true;
    while (solved && step.substep < step.substeps) {
        double w = frequency.angular = 2 * M_PI * frequency.current;

        eslog::info(" ============================================================================================= \n");
        eslog::info(" = LOAD STEP %2d                                                         FREQUENCY %10.0f = \n", step.loadstep + 1, frequency.current);
        eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
        double start = eslog::time();
        assembler.evaluate(step, frequency, K, M, C, re.f, im.f, nullptr, nullptr, re.dirichlet, im.dirichlet);
        eslog::checkpointln("SIMULATION: PHYSICS ASSEMBLED");
        storeSystem(step);

        solver->A->set(0);
        solver->A->copy(K, s1, s1)->add(-w * w, M, s1, s1);
        solver->A->copy(K, s2, s2)->add(-w * w, M, s2, s2);
        solver->A->updated = true;
        solver->b->copy(re.f, s1);
        solver->b->copy(im.f, s2);
        solver->b->updated = true;
        solver->dirichlet->copy(re.dirichlet, s1);
        solver->dirichlet->copy(im.dirichlet, s2);
        solver->dirichlet->updated = true;
        eslog::info("       = ----------------------------------------------------------------------------- = \n");
        eslog::info("       = SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

        solver->update(step);
        eslog::checkpointln("SIMULATION: LINEAR SYSTEM UPDATED");
        solved = solver->solve(step);
        eslog::checkpointln("SIMULATION: LINEAR SYSTEM SOLVED");

        double solution = eslog::time();
        re.x->copy(solver->x, s1);
        im.x->copy(solver->x, s2);
        storeSolution(step);
        assembler.updateSolution(step, re.x, im.x);
        info::mesh->output->updateSolution(step, frequency);
        eslog::info("       = PROCESS SOLUTION                                                   %8.3f s = \n", eslog::time() - solution);
        eslog::info("       = ----------------------------------------------------------------------------- = \n");

        eslog::info(" ====================================================================== solved in %8.3f s = \n\n", eslog::time() - start);
        eslog::checkpointln("SIMULATION: SOLUTION PROCESSED");
        frequency.current += frequency.shift;
        ++step.substep;
    }
    return solved;
}

void StructuralMechanicsHarmonicRealLinear::storeSystem(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{K, M, C, f.re, f.im, dirichlet.re}\n");
        K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
        M->store(utils::filename(utils::debugDirectory(step) + "/scheme", "M").c_str());
        C->store(utils::filename(utils::debugDirectory(step) + "/scheme", "C").c_str());
        re.f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f.re").c_str());
        im.f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f.im").c_str());
        re.dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet.re").c_str());
        im.dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet.im").c_str());
    }
}

void StructuralMechanicsHarmonicRealLinear::storeSolution(step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: scheme/{x.re, x.im}\n");
        re.x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x.re").c_str());
        im.x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x.im").c_str());
    }
}
