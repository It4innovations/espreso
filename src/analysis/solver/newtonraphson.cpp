
#include "newtonraphson.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/physicssolver/nonlinear.h"

#include "analysis/assembler/module/heattransfer.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"

using namespace espreso;

AX_NewtonRaphson::AX_NewtonRaphson(NonLinearSolverConfiguration &configuration)
: configuration(configuration), U(nullptr), R(nullptr)
{

}

AX_NewtonRaphson::~AX_NewtonRaphson()
{

}

void AX_NewtonRaphson::init(AX_LinearSystem<double> *system)
{
	eslog::info(" == NON-LINEAR SOLVER                                                        NEWTON-RAPHSON == \n");

	U = system->solver.x->copyPattern();
	R = system->solver.x->copyPattern();
	system->solver.A->commit();
}

bool AX_NewtonRaphson::checkTemp(step::Step &step, AX_HeatTransfer &assembler, AX_SteadyState &scheme, AX_LinearSystem<double> *system)
{
	double solutionNumerator = system->solver.x->norm();
	system->solver.x->add(1, U);
	scheme.extractSolution(step, system);

	double solutionDenominator = std::max(system->solver.x->norm(), 1e-3);
	double norm = solutionNumerator / solutionDenominator;

	if (norm > configuration.requested_first_residual) {
		eslog::info("     == TEMPERATURE NORM, CRITERIA                            %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.requested_first_residual);
	} else {
		eslog::info("     == TEMPERATURE NORM, CRITERIA [CONVERGED]                %.5e / %.5e == \n", solutionNumerator, solutionDenominator * configuration.requested_first_residual);
	}

	assembler.updateSolution();
	return !(norm > configuration.requested_first_residual);
}

bool AX_NewtonRaphson::run(step::Step &step, step::Time &time, AX_HeatTransfer &assembler, AX_SteadyState &scheme, AX_LinearSystem<double> *system)
{
	eslog::info("     ===================================================================================== \n");
	eslog::info("     ==  NEWTON RAPHSON CONVERGENCE CRITERIA                                            == \n");
	if (configuration.check_first_residual) {
		eslog::info("     ==  - TEMPERATURE RESIDUAL                                                    TRUE == \n");
	} else {
		eslog::info("     ==  - TEMPERATURE RESIDUAL                                                   FALSE == \n");
	}
	if (configuration.check_second_residual) {
		eslog::info("     ==  - HEAT RESIDUAL                                                           TRUE == \n");
	} else {
		eslog::info("     ==  - HEAT RESIDUAL                                                          FALSE == \n");
	}
	eslog::info("     == ------------------------------------------------------------------------------- == \n");

	system->info();
	eslog::info("     == ------------------------------------------------------------------------------- == \n");
	eslog::info("     == INITIAL STEP                                                                    ==     \n");

	step.iteration = 0;
	assembler.evaluate();
	scheme.composeSystem(step, system);
//	assembler.fillDirichlet(*system->solver.dirichlet);

	system->update(step, assembler);
	system->solve(step);

	scheme.extractSolution(step, system);
	assembler.updateSolution();

	// iterations
	while (step.iteration++ < configuration.max_iterations) {
		eslog::info("     == ------------------------------------------------------------------------------- == \n");
		eslog::info("\n     == EQUILIBRIUM ITERATION                                                      [%2d] == \n", step.iteration);

		U->fillData(system->solver.x);

		assembler.evaluate();
		scheme.composeSystem(step, system);

		system->solver.A->apply(1, system->solver.x, 0, R);
		system->solver.b->add(-1, R);
		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: scheme/{R}\n");
			R->store(utils::filename(utils::debugDirectory(step) + "/scheme", "R").c_str());
		}

//		assembler.fillDirichlet(*system->solver.dirichlet);
//		U->addTo(-1, system->solver.dirichlet);
		system->solver.dirichlet->add(-1, U);

		system->update(step, assembler);
		system->solve(step);

		if (checkTemp(step, assembler, scheme, system)) {
			break;
		}
	}

	return true;
}
