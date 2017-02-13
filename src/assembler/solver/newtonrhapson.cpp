
#include "newtonrhapson.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"
#include "../physics/physics.h"

#include "../../configuration/physics/nonlinearsolver.h"

#include "../../solver/generic/LinearSolver.h"

using namespace espreso;

NewtonRhapson::NewtonRhapson(
		Mesh *mesh,
		std::vector<Physics*> &physics,
		std::vector<Instance*> &instances,
		std::vector<LinearSolver*> &linearSolvers,
		store::ResultStore* store,
		const NonLinearSolver &configuration)
: Solver(mesh, physics, instances, linearSolvers, store), _configuration(configuration)
{

}

void NewtonRhapson::run(Step &step)
{
	if (!_configuration.convergence_parameters.temperature && _configuration.convergence_parameters.heat) {
		ESINFO(GLOBAL_ERROR) << "It is not possible to turn off the both 'temperature' and 'heat' convergence.";
	}
	step.solver = 0;

	assembleStiffnessMatrices(step);
	assembleB1(step);
	makeStiffnessMatricesRegular(step);
	assembleB0(step);

	initLinearSolver();
	startLinearSolver();
	storeSolution(step);

	double temperatureResidual = _configuration.convergence_parameters.temperature_residual;
	double heatResidual = _configuration.convergence_parameters.heat_residual;
	if (_configuration.convergence_parameters.temperature) {
		temperatureResidual *= 10;
	}
	if (_configuration.convergence_parameters.heat) {
		heatResidual *= 10;
	}

	while (
		step.solver++ < _configuration.max_iterations &&
		(temperatureResidual > _configuration.convergence_parameters.temperature_residual ||
		heatResidual > _configuration.convergence_parameters.heat_residual)) {

		std::vector<std::vector<double> > T = instances.front()->primalSolution;

		instances.front() = new Instance(instances.front()->domains);
		instances.front()->DOFs = physics.front()->instance()->DOFs;
		instances.front()->primalSolution = physics.front()->instance()->primalSolution;
		instances.front()->solutions = physics.front()->instance()->solutions;
		physics.front()->instance()->solutions.resize(0);
		linearSolvers.front() = new LinearSolver(linearSolvers.front()->configuration, linearSolvers.front()->physics, linearSolvers.front()->constraints);
		physics.front()->_instance = instances.front();

		assembleStiffnessMatrices(step);
		subtractResidualForces(step);
		assembleB1(step);

		subtractSolutionFromB1c(step);

		makeStiffnessMatricesRegular(step);

		initLinearSolver();
		startLinearSolver();

		if (_configuration.convergence_parameters.temperature) {
			temperatureResidual = physics.front()->computeNormOfSolution();
		}
		addToSolution(physics.front(), T);
		if (_configuration.convergence_parameters.temperature) {
			temperatureResidual /= physics.front()->computeNormOfSolution();
		}
		storeSolution(step);
	}

	finalizeLinearSolver();
}



