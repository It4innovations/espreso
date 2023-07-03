

#include <analysis/physics/physics.h>
#include "structuralmechanics.steadystate.nonlinear.h"

#include "analysis/linearsystem/feti/fetisystem.h"
#include "analysis/linearsystem/direct/mklpdsssystem.h"
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
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, U{}, R{}, f{}, x{}, dirichlet{}, system{}
{

}

StructuralMechanicsSteadyStateNonLinear::~StructuralMechanicsSteadyStateNonLinear()
{
	if (system) { delete system; }
	if (K) { delete K; }
	if (U) { delete U; }
	if (R) { delete R; }
	if (f) { delete f; }
	if (x) { delete x; }
	if (dirichlet) { delete dirichlet; }
}

void StructuralMechanicsSteadyStateNonLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                        PLASTICITY STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                            STRUCTURAL MECHANICS == \n");
	eslog::info(" ============================================================================================= \n");

	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::TIME);
	if (configuration.nonlinear_solver.stepping == NonLinearSolverConfiguration::STEPPINGG::FALSE) {
		configuration.nonlinear_solver.substeps = 1;
	}
}

void StructuralMechanicsSteadyStateNonLinear::run(step::Step &step)
{
	switch (configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new FETISystem<StructuralMechanicsSteadyStateNonLinear>(this); break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new MKLPDSSSystem<StructuralMechanicsSteadyStateNonLinear>(this); break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    break;
	}
	system->setMapping(K = system->assembler.A->copyPattern());
	system->setMapping(R = system->solver.x->copyPattern());
	system->setMapping(f = system->assembler.b->copyPattern());
	system->setMapping(x = system->assembler.x->copyPattern());
	system->setDirichletMapping(dirichlet = system->assembler.dirichlet->copyPattern());
	U = system->solver.x->copyPattern();
	system->solver.A->commit();

	assembler.connect(K, nullptr, nullptr, f, R, dirichlet);

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
	system->set(step);
	eslog::info(" ============================================================================================= \n\n");
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM SET");

	eslog::info(" ============================================================================================= \n");
	eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step::step.loadstep + 1, time.current);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	eslog::info("      =================================================================================== \n");
	eslog::info("      ==  NEWTON RAPHSON CONVERGENCE CRITERIA                                          == \n");
	if (configuration.nonlinear_solver.check_first_residual) {
		eslog::info("      ==  - DISPLACEMENT RESIDUAL                                                 TRUE == \n");
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
	for (step.substep = 0; step.substep < configuration.nonlinear_solver.substeps; ++step.substep) {
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

		assembler.evaluate(step, time, K, nullptr, nullptr, f, nullptr, dirichlet);
		storeSystem(step);
		system->solver.A->copy(K);
		system->solver.b->copy(f);
		system->solver.dirichlet->copy(dirichlet);
		eslog::info("      == ----------------------------------------------------------------------------- == \n");
		eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

		system->update(step);
		system->solve(step);

		double solution = eslog::time();
		x->copy(system->solver.x);
		storeSolution(step);
		assembler.updateSolution(x);
		eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
		eslog::info("      == ----------------------------------------------------------------------------- == \n");

		// iterations
		while (step.iteration++ < configuration.nonlinear_solver.max_iterations) {
			eslog::info("\n      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

			start = eslog::time();
			U->copy(system->solver.x);
			assembler.evaluate(step, time, K, nullptr, nullptr, f, R, dirichlet);
			storeSystem(step);
			system->solver.A->copy(K);
			system->solver.b->copy(f);
			system->solver.b->add(-1, R);
			system->solver.dirichlet->copy(dirichlet);
			system->solver.dirichlet->add(-1, U);

			eslog::info("      == ----------------------------------------------------------------------------- == \n");
			eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s == \n", eslog::time() - start);

			system->update(step);
			system->solve(step);

			if (checkDisplacement(step)) {
				break;
			}
		}
		info::mesh->output->updateSolution(step, time);
	}
}

bool StructuralMechanicsSteadyStateNonLinear::checkDisplacement(step::Step &step)
{
	double solution = eslog::time();
	double solutionNumerator = system->solver.x->norm();
	system->solver.x->add(1, U);
	x->copy(system->solver.x);

	double solutionDenominator = std::max(system->solver.x->norm(), 1e-3);
	double norm = solutionNumerator / solutionDenominator;

	eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");

	if (norm > configuration.nonlinear_solver.requested_first_residual) {
		eslog::info("      == DISPLACEMENT NORM, CRITERIA                         %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.nonlinear_solver.requested_first_residual);
		assembler.nextIteration(x);
		return false;
	} else {
		eslog::info("      == DISPLACEMENT NORM, CRITERIA [CONVERGED]             %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.nonlinear_solver.requested_first_residual);
		eslog::info("      =================================================================================== \n\n");
		assembler.updateSolution(x);
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
		dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet").c_str());
	}
}

void StructuralMechanicsSteadyStateNonLinear::storeSolution(step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{x}\n");
		x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
	}
}
