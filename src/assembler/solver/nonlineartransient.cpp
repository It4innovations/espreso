
#include "nonlineartransient.h"


#include "../../configuration/physics/transientsolver.h"
#include "../../configuration/physics/nonlinearsolver.h"

#include "../step.h"
#include "../instance.h"
#include "../solution.h"
#include "../../mesh/structures/elementtypes.h"
#include "../physics/physics.h"
#include "../../basis/logging/logging.h"
#include "../../basis/logging/constants.h"

#include <cmath>

using namespace espreso;

size_t NonlinearTransient::offset = -1;
size_t NonlinearTransient::lastStep = -1;

NonlinearTransient::NonlinearTransient(
		Mesh *mesh,
		Physics* physics,
		FETISolver* linearSolver,
		output::Store* store,
		const TransientSolver &configuration,
		double duration)
: Solver("TRANSIENT FIRST ORDER IMPLICIT", mesh, physics, linearSolver, store, duration, Matrices::NONE),
  _configuration(configuration)
{

}

void NonlinearTransient::run(Step &step)
{
	ESINFO(PROGRESS1) << "Run " << _name << " solver for " << physics->name();

	double alpha = 0;
	switch (_configuration.method) {
	case TransientSolver::METHOD::CRANK_NICOLSON:
		alpha = 0.5;
		break;
	case TransientSolver::METHOD::GALERKIN:
		alpha = 2 / 3;
		break;
	case TransientSolver::METHOD::BACKWARD_DIFF:
		alpha = 1;
		break;
	case TransientSolver::METHOD::USER:
		alpha = _configuration.alpha;
		if (alpha <= 0 || alpha > 1) {
			ESINFO(GLOBAL_ERROR) << "Alpha has to be from interval (0, 1>.";
		}
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not supported first order implicit solver method.";
	}

	preprocessData(step);

	if (offset == (size_t)-1) {
		lastStep = step.step;
		offset = instance->solutions.size();
		instance->solutions.push_back(new Solution(*_mesh, "trans_U" , ElementType::NODES, physics->pointDOFs()));
		instance->solutions.push_back(new Solution(*_mesh, "trans_dU", ElementType::NODES, physics->pointDOFs()));
		instance->solutions.push_back(new Solution(*_mesh, "trans_V" , ElementType::NODES, physics->pointDOFs()));
		instance->solutions.push_back(new Solution(*_mesh, "trans_X" , ElementType::NODES, physics->pointDOFs()));
		instance->solutions.push_back(new Solution(*_mesh, "trans_Y" , ElementType::NODES, physics->pointDOFs()));
	}

	std::vector<std::vector<double> > &u     = instance->solutions[offset + SolutionIndex::SOLUTION]->innerData();
	std::vector<std::vector<double> > deltaU = instance->solutions[offset + SolutionIndex::DELTA]->innerData();
	std::vector<std::vector<double> > v      = instance->solutions[offset + SolutionIndex::VELOCITY]->innerData();
	std::vector<std::vector<double> > x      = instance->solutions[offset + SolutionIndex::X]->innerData();
	std::vector<std::vector<double> > y      = instance->solutions[offset + SolutionIndex::Y]->innerData();

	u = instance->primalSolution;
	if (lastStep + 1 != step.step) {
		lastStep = step.step;
		for (size_t i = 0; i < v.size(); i++) {
			std::fill(v[i].begin(), v[i].end(), 0);
		}
	}

	double startTime = step.currentTime;
	step.timeStep = _configuration.time_step;
	step.iteration = 0;
	step.currentTime += step.timeStep;
	bool timeDependent = physics->isMatrixTimeDependent(step);

	for (step.substep = 0; step.currentTime <= startTime + _duration + 1e-8; step.substep++) {
		ESINFO(PROGRESS2) << _name << " iteration " << step.substep + 1 << "[" << step.currentTime << "s]";
		setEmptyRegularization(step, Matrices::K);
		if (step.substep == 0 || timeDependent) {
			updateMatrices(step, Matrices::K | Matrices::M | Matrices::f);
			composeGluing(step, Matrices::B1 | Matrices::B0); // TODO: B0 without kernels
			sum(instance->K, 1 / (alpha * _configuration.time_step), instance->M, "K += (1 / " + ASCII::alpha + ASCII::DELTA + "t) * M");
		} else {
			updateMatrices(step, Matrices::f);
			composeGluing(step, Matrices::B1); // TODO: B0 without kernels
		}

		sum(x, 1 / (alpha * _configuration.time_step), u, (1 - alpha) / alpha, v, "x = (1 / " + ASCII::alpha + ASCII::DELTA + ") * u + (1 - " + ASCII::alpha + ") / " + ASCII::alpha + " * v");

		multiply(y, instance->M, x, "y = M * x");
		sum(instance->f, 1, instance->f, 1, y, "f += y");

		if (step.substep) {
			updateLinearSolver(step, Matrices::f | Matrices::B1c);
		} else {
			initLinearSolver(step);
		}
		runLinearSolver(step);

		/////////////////////////////////////////////////



		std::vector<std::vector<double> > T;
		std::vector<std::vector<double> > F_ext;
		std::vector<std::vector<double> > f_R_BtLambda;

		int cumiter = 0;

		double temperatureResidual = 10 * 1e-3; //_configuration.convergenceParameters().requestedSolution();
		double temperatureResidual_first = 0;
		double temperatureResidual_second = 0;

		double heatResidual;
		double heatResidual_first = 0;
		double heatResidual_second = 0;

		while (step.iteration++ < 100) { //_configuration.max_iterations) {

			cumiter +=1;
			T = physics->instance()->primalSolution;

			step.timeIntegrationConstantM = 1 / (alpha * _configuration.time_step);
			updateMatrices(step, Matrices::M | Matrices::K | Matrices::f | Matrices::R, physics->instance()->solutions);
			sum(instance->K, 1 / (alpha * _configuration.time_step), instance->M, "K += (1 / " + ASCII::alpha + ASCII::DELTA + "t) * M");

			sum(x, 1 / (alpha * _configuration.time_step), u, (1 - alpha) / alpha, v, "x = (1 / " + ASCII::alpha + ASCII::DELTA + ") * u + (1 - " + ASCII::alpha + ") / " + ASCII::alpha + " * v");

			multiply(y, instance->M, x, "y = M * x");
			sum(instance->f, 1, instance->f, 1, y, "f += y");

			sum(instance->f, 1, instance->f, -1, instance->R, "f = f - R");

			composeGluing(step, Matrices::B1);
			subtractDirichlet();

			updateLinearSolver(step, Matrices::K | Matrices::f | Matrices::B1c);

			runLinearSolver(step);

			ESINFO(CONVERGENCE) <<  "    LINEAR_SOLVER_OUTPUT: SOLVER = " << "PCG" <<   " N_ITERATIONS = " << "1" << "  " ;

			temperatureResidual_first = sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));

			sum(instance->primalSolution, 1, instance->primalSolution, 1, T, "u = " + ASCII::DELTA + "u + u");
			temperatureResidual_second = sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));
			if (temperatureResidual_second < 1e-3) {
				temperatureResidual_second = 1e-3;
			}
			temperatureResidual = temperatureResidual_first / temperatureResidual_second;

			if ( temperatureResidual > 1e-3) { //_configuration.convergenceParameters().requestedSolution()) {
				ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  temperatureResidual_first << "  CRITERION_VALUE = " << temperatureResidual_second * 1e-3; //_configuration.convergenceParameters().requestedSolution() ;
			} else {
				ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  temperatureResidual_first << "  CRITERION_VALUE = " << temperatureResidual_second * 1e-3 /*_configuration.convergenceParameters().requestedSolution()*/ <<  " <<< CONVERGED >>>" ;
				break;
			}

			storeSubSolution(step);
			processSolution(step);
		}
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << step.iteration + 1 ;
		ESINFO(CONVERGENCE) <<  " >> SUBSTEP " << step.substep + 1 << " IS DONE /" << " CUMULATIVE ITERATION NUMBER = " << cumiter + 1 ;
		ESINFO(CONVERGENCE) <<  " ";
		step.iteration = 0;

		/////////////////////////////////////////////////


		sum(deltaU, 1, instance->primalSolution, -1, u, ASCII::DELTA + "u = u_i - u_i_1");
		sum(v, 1 / (alpha * _configuration.time_step), deltaU, - (1 - alpha) / alpha, v, "v = (1 / " + ASCII::alpha + ASCII::DELTA + "t) * " + ASCII::DELTA + "u - (1 - " + ASCII::alpha + ") / " + ASCII::alpha + " * v");
		u = instance->primalSolution;

		processSolution(step);
		storeSolution(step);
		step.currentTime += step.timeStep;
	}
	finalizeLinearSolver(step);
}

void NonlinearTransient::init(Step &step)
{

}

void NonlinearTransient::preprocess(Step &step)
{
}

void NonlinearTransient::solve(Step &step)
{

}

void NonlinearTransient::postprocess(Step &step)
{
}

void NonlinearTransient::finalize(Step &step)
{
}



