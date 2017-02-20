
#include "newtonrhapson.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"
#include "../physics/physics.h"

#include "../../basis/utilities/utils.h"
#include "../../configuration/physics/nonlinearsolver.h"
#include "../../solver/generic/LinearSolver.h"

using namespace espreso;

NewtonRhapson::NewtonRhapson(
		Mesh *mesh,
		Physics* physics,
		LinearSolver* linearSolver,
		store::ResultStore* store,
		const NonLinearSolverBase &configuration)
: Solver("NEWTON RHAPSON", mesh, physics, linearSolver, store), _configuration(configuration)
{

}

void NewtonRhapson::run(Step &step)
{
	ESINFO(PROGRESS1) << "Run " << _name << " solver for " << physics->name();

	if (!_configuration.convergenceParameters().checkSolution() && _configuration.convergenceParameters().checkResidual()) {
		ESINFO(GLOBAL_ERROR) << "It is not possible to turn off the both 'temperature' and 'heat' convergence.";
	}
	step.solver = 0;

	assembleMatrices(step, Matrices::K | Matrices::f);
	composeGluing(step, Matrices::B1);
	regularizeMatrices(step, Matrices::K);
	composeGluing(step, Matrices::B0);

	initLinearSolver();
	runLinearSolver();
	processSolution(step);

	double temperatureResidual = _configuration.convergenceParameters().requestedSolution();
	double heatResidual = _configuration.convergenceParameters().requestedResidual();
	if (_configuration.convergenceParameters().checkSolution()) {
		temperatureResidual *= 10;
	}
	if (_configuration.convergenceParameters().checkResidual()) {
		heatResidual *= 10;
	}

	std::vector<std::vector<double> > T;
	std::vector<std::vector<double> > F_ext;
	while (
		step.solver++ < _configuration.max_iterations &&
		(temperatureResidual > _configuration.convergenceParameters().requestedSolution() ||
		heatResidual > _configuration.convergenceParameters().requestedResidual())) {

		T = physics->instance()->primalSolution;

		updateMatrices(step, Matrices::K | Matrices::f | Matrices::R, physics->instance()->solutions);

		if (_configuration.line_search) {
			F_ext = physics->instance()->f;
		}
		if (_configuration.convergenceParameters().checkResidual()) {
			heatResidual = physics->sumSquares(physics->instance()->f, Physics::SumOperation::SUM);
		}
		updateVector(step, Matrices::f, Matrices::R, 1, -1);
		if (_configuration.convergenceParameters().checkResidual()) {
			heatResidual += physics->sumSquares(physics->instance()->f, Physics::SumOperation::SUM, Physics::SumRestriction::DIRICHLET, step.load);
		}

		composeGluing(step, Matrices::B1);
		updateVector(step, Matrices::B1c, Matrices::primar, 1, -1);
		regularizeMatrices(step, Matrices::K);
		updateLinearSolver(Matrices::K | Matrices::f | Matrices::B1c);
		runLinearSolver();

		if (_configuration.line_search) {
			lineSearch(T, physics->instance()->primalSolution, F_ext, physics, step);
		}
		if (_configuration.convergenceParameters().checkSolution()) {
			temperatureResidual = sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));
		}
		sumVectors(physics->instance()->primalSolution, T, physics->instance()->primalSolution);
		updateVector(step, Matrices::primar, T, 1, 1);
		if (_configuration.convergenceParameters().checkSolution()) {
			temperatureResidual /= sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));;
		}

		processSolution(step);

		if (_configuration.convergenceParameters().checkResidual()) {
			assembleMatrices(step, Matrices::K | Matrices::R | Matrices::f);
			// double R;
			// |R| = f - R
			// R = physics.front()->sumSquares(physics.front()->instance()->f, Physics::SumOperation::SUM);
			// |R| -= BtLambda
			// TODO
		}
	}
	storeSolution(step);
	finalizeLinearSolver();
}



