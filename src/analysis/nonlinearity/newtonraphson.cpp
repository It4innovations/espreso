
#include "newtonraphson.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/physicssolver/nonlinear.h"

#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/structuralmechanics.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"

using namespace espreso;

NewtonRaphson::NewtonRaphson(NonLinearSolverConfiguration &configuration)
: configuration(configuration), U(nullptr), R(nullptr)
{

}

NewtonRaphson::~NewtonRaphson()
{

}

void NewtonRaphson::init(LinearSystem<double> *system)
{
	U = system->solver.x->copyPattern();
	R = system->solver.x->copyPattern();
	system->solver.A->commit();
}

bool NewtonRaphson::run(step::Step &step, step::Time &time, HeatTransfer &assembler, SteadyState &scheme, LinearSystem<double> *system)
{
	eslog::info("      =================================================================================== \n");
	eslog::info("      ==  NEWTON RAPHSON CONVERGENCE CRITERIA                                          == \n");
	if (configuration.check_first_residual) {
		eslog::info("      ==  - TEMPERATURE RESIDUAL                                                  TRUE == \n");
	} else {
		eslog::info("      ==  - TEMPERATURE RESIDUAL                                                 FALSE == \n");
	}
	if (configuration.check_second_residual) {
		eslog::info("      ==  - HEAT RESIDUAL                                                         TRUE == \n");
	} else {
		eslog::info("      ==  - HEAT RESIDUAL                                                        FALSE == \n");
	}
	eslog::info("      =================================================================================== \n");

	eslog::info("      ==                                                                  INITIAL STEP ==     \n");

	double start = eslog::time();
	step.iteration = 0;
	assembler.evaluate(scheme, time);
	scheme.composeSystem(step, system);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");
	eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

	system->update(step);
	system->solve(step);

	double solution = eslog::time();
	scheme.extractSolution(step, system);
	assembler.updateSolution(scheme);
	eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");

	// iterations
	while (step.iteration++ < configuration.max_iterations) {
		eslog::info("\n      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

		start = eslog::time();
		U->copy(system->solver.x);
		assembler.evaluate(scheme, time);
		scheme.composeSystem(step, system);

		system->solver.A->apply(1, system->solver.x, 0, R);
		system->solver.b->add(-1, R);

		// why not to set dirichlet to 0 for all iterations??
		system->solver.dirichlet->add(-1, U);

		eslog::info("      == ----------------------------------------------------------------------------- == \n");
		eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s == \n", eslog::time() - start);

		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: scheme/{R}\n");
			R->store(utils::filename(utils::debugDirectory(step) + "/scheme", "R").c_str());
		}

		system->update(step);
		system->solve(step);

		if (checkTemp(step, assembler, scheme, system)) {
			break;
		}
	}

	return true;
}

bool NewtonRaphson::checkTemp(step::Step &step, HeatTransfer &assembler, SteadyState &scheme, LinearSystem<double> *system)
{
	double solution = eslog::time();
	double solutionNumerator = system->solver.x->norm();
	system->solver.x->add(1, U);
	scheme.extractSolution(step, system);

	double solutionDenominator = std::max(system->solver.x->norm(), 1e-3);
	double norm = solutionNumerator / solutionDenominator;

	eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");

	if (norm > configuration.requested_first_residual) {
		eslog::info("      == TEMPERATURE NORM, CRITERIA                          %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.requested_first_residual);
	} else {
		eslog::info("      == TEMPERATURE NORM, CRITERIA [CONVERGED]              %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.requested_first_residual);
	}

	assembler.updateSolution(scheme);
	return !(norm > configuration.requested_first_residual);
}

bool NewtonRaphson::run(step::Step &step, step::Time &time, StructuralMechanics &assembler, SteadyState &scheme, LinearSystem<double> *system)
{
	eslog::info("      =================================================================================== \n");
	eslog::info("      ==  NEWTON RAPHSON CONVERGENCE CRITERIA                                          == \n");
	if (configuration.check_first_residual) {
		eslog::info("      ==  - DISPLACEMENT RESIDUAL                                                 TRUE == \n");
	} else {
		eslog::info("      ==  - DISPLACEMENT RESIDUAL                                                FALSE == \n");
	}
	if (configuration.check_second_residual) {
		eslog::info("      ==  - STRESS RESIDUAL                                                       TRUE == \n");
	} else {
		eslog::info("      ==  - STRESS RESIDUAL                                                      FALSE == \n");
	}
	eslog::info("      =================================================================================== \n");

	eslog::info("      ==                                                                  INITIAL STEP ==     \n");

	double start = eslog::time();
	step.iteration = 0;
	assembler.evaluate(scheme, time);
	scheme.composeSystem(step, system);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");
	eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s = \n", eslog::time() - start);

	system->update(step);
	system->solve(step);

	double solution = eslog::time();
	scheme.extractSolution(step, system);
	assembler.updateSolution(scheme);
	eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");

	// iterations
	while (step.iteration++ < configuration.max_iterations) {
		eslog::info("\n      ==                                                    %3d. EQUILIBRIUM ITERATION == \n", step.iteration);

		start = eslog::time();
		U->copy(system->solver.x);
		assembler.evaluate(scheme, time);
		scheme.composeSystem(step, system);

		system->solver.A->apply(1, system->solver.x, 0, R);
		system->solver.b->add(-1, R);

		// why not to set dirichlet to 0 for all iterations??
		system->solver.dirichlet->add(-1, U);

		eslog::info("      == ----------------------------------------------------------------------------- == \n");
		eslog::info("      == SYSTEM ASSEMBLY                                                    %8.3f s == \n", eslog::time() - start);

		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: scheme/{R}\n");
			R->store(utils::filename(utils::debugDirectory(step) + "/scheme", "R").c_str());
		}

		system->update(step);
		system->solve(step);

		if (checkDisplacement(step, assembler, scheme, system)) {
			break;
		}
	}

	return true;
}

bool NewtonRaphson::checkDisplacement(step::Step &step, StructuralMechanics &assembler, SteadyState &scheme, LinearSystem<double> *system)
{
	double solution = eslog::time();
	double solutionNumerator = system->solver.x->norm();
	system->solver.x->add(1, U);
	scheme.extractSolution(step, system);

	double solutionDenominator = std::max(system->solver.x->norm(), 1e-3);
	double norm = solutionNumerator / solutionDenominator;

	eslog::info("      == PROCESS SOLUTION                                                   %8.3f s == \n", eslog::time() - solution);
	eslog::info("      == ----------------------------------------------------------------------------- == \n");

	if (norm > configuration.requested_first_residual) {
		eslog::info("      == DISPLACEMENT NORM, CRITERIA                         %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.requested_first_residual);
	} else {
		eslog::info("      == DISPLACEMENT NORM, CRITERIA [CONVERGED]             %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.requested_first_residual);
	}

	assembler.updateSolution(scheme);
	return !(norm > configuration.requested_first_residual);
}

