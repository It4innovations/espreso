
#include "physics.h"
#include "heat.steadystate.nonlinear.h"

#include "analysis/builder/uniformbuilder.direct.h"
#include "analysis/builder/uniformbuilder.feti.h"
#include "analysis/linearsystem/mklpdsssolver.h"
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
: settings(settings), configuration(configuration), assembler{nullptr, settings, configuration}, K{}, U{}, R{}, f{}, x{}, dirichlet{}, builder{}, solver{}

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
	if (builder) { delete builder; }
	if (solver) { delete solver; }
}

void HeatSteadyStateNonLinear::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                        NON-LINEAR STEADY STATE == \n");
	eslog::info(" == PHYSICS                                                                   HEAT TRANSFER == \n");
	eslog::info(" ============================================================================================= \n");

	assembler.analyze();
	info::mesh->output->updateMonitors(step::TYPE::TIME);

	Matrix_Shape shape = Matrix_Shape::UPPER;
	Matrix_Type type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	for (auto mat = settings.material_set.begin(); mat != settings.material_set.end(); ++mat) {
		if (settings.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
			shape = Matrix_Shape::FULL;
			Matrix_Type type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
		}
	}
	if (configuration.translation_motions.size()) {
		shape = Matrix_Shape::FULL;
		type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
	}

	switch (configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:
		builder = new UniformBuilderFETI<double>(configuration.temperature, 1, shape);
		break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS:
		builder = new UniformBuilderDirect<double>(configuration.temperature, 1, shape);
		solver = new MKLPDSSLinearSystemSolver<double>(configuration.mklpdss);
		break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    break;
	}

	builder->fillMatrix(solver->A, type, shape);
	builder->fillVector(solver->b);
	builder->fillVector(solver->x);
	builder->fillDirichlet(solver->dirichlet);

	K = solver->A->copyPattern();
	R = solver->b->copyPattern();
	U = solver->b->copyPattern();
	f = solver->b->copyPattern();
	x = solver->x->copyPattern();
	dirichlet = solver->dirichlet->copyPattern();

	builder->fillMatrixMap(K);
	builder->fillVectorMap(R);
	builder->fillVectorMap(f);
	builder->fillDirichletMap(dirichlet);
	eslog::checkpointln("SIMULATION: LINEAR SYSTEM BUILT");
}

void HeatSteadyStateNonLinear::run(step::Step &step)
{
	time.shift = configuration.duration_time;
	time.start = 0;
	time.current = configuration.duration_time;
	time.final = configuration.duration_time;

	assembler.connect(K, nullptr, f, R, dirichlet);

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
	eslog::info(" = LOAD STEP %2d                                                              TIME %10.4f = \n", step::step.loadstep + 1, time.current);
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
	assembler.evaluate(time, K, nullptr, f, nullptr, dirichlet);
	storeSystem(step);
	solver->A->copy(K);
	solver->b->copy(f);
	solver->dirichlet->copy(dirichlet);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");
	eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

	solver->update(step);
	solver->solve(step);

	double solution = eslog::time();
	x->copy(solver->x);
	storeSolution(step);
	assembler.updateSolution(x);
	eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");

	// iterations
	while (step.iteration++ < configuration.nonlinear_solver.max_iterations) {
		eslog::info("\n      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

		start = eslog::time();
		U->copy(solver->x);
		assembler.evaluate(time, K, nullptr, f, R, dirichlet);
		storeSystem(step);
		solver->A->copy(K);
		solver->b->copy(f);
		solver->b->add(-1, R);
		solver->dirichlet->copy(dirichlet);
		solver->dirichlet->add(-1, U);

		eslog::info("      == ----------------------------------------------------------------------------- == \n");
		eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s == \n", eslog::time() - start);

		solver->update(step);
		solver->solve(step);

		if (checkTemp(step)) {
			break;
		}
	}
	info::mesh->output->updateSolution(step, time);
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

	assembler.updateSolution(x);
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


