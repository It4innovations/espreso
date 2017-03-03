
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
		const NonLinearSolverBase &configuration,
		Matrices restriction)
: Solver("NEWTON RHAPSON", mesh, physics, linearSolver, store, restriction), _configuration(configuration)
{

}

void NewtonRhapson::run(Step &step)
{
	ESINFO(PROGRESS1) << "Run " << _name << " solver for " << physics->name();

	init(step);
	preprocess(step);
	solve(step);
	postprocess(step);
	finalize(step);
}

void NewtonRhapson::init(Step &step)
{
	assembleMatrices(step, Matrices::K | Matrices::f);
	size_t substeps = _configuration.stepping == NonLinearSolverBase::STEPPINGG::TRUE ? _configuration.substeps : 1;
	multiply(step, Matrices::f, (step.iteration + 1) / (double)substeps);
	composeGluing(step, Matrices::B1);
	multiply(step, Matrices::B1c, (step.iteration + 1) / (double)substeps);
	regularizeMatrices(step, Matrices::K);
	composeGluing(step, Matrices::B0);
}

void NewtonRhapson::preprocess(Step &step)
{
	initLinearSolver();
}

void NewtonRhapson::solve(Step &step)
{
	if (!_configuration.convergenceParameters().checkSolution() && !_configuration.convergenceParameters().checkResidual()) {
		ESINFO(GLOBAL_ERROR) << "It is not possible to turn off the both 'temperature' and 'heat' convergence.";
	}

	runLinearSolver();
	processSolution(step);
	storeSubSolution(step);

	std::vector<std::vector<double> > T;
	std::vector<std::vector<double> > F_ext;
	std::vector<std::vector<double> > f_R_BtLambda;


	size_t substeps = _configuration.stepping == NonLinearSolverBase::STEPPINGG::TRUE ? _configuration.substeps : 1;
	for (step.iteration = 0; step.iteration < substeps; step.iteration++) {

		step.substep = 0;
		double temperatureResidual = 10 * _configuration.convergenceParameters().requestedSolution();
		double heatResidual;

		while (
			step.substep++ < _configuration.max_iterations &&
			temperatureResidual > _configuration.convergenceParameters().requestedSolution()) {

			T = physics->instance()->primalSolution;

			if (_configuration.method == NonLinearSolverBase::METHOD::MODIFIED_NEWTON_RHAPSON && step.substep > 1) {
				updateMatrices(step, Matrices::f | Matrices::R, physics->instance()->solutions);
			} else {
				updateMatrices(step, Matrices::K | Matrices::f | Matrices::R, physics->instance()->solutions);
			}
			multiply(step, Matrices::f, (step.iteration + 1) / (double)substeps);

			if (_configuration.line_search) {
				F_ext = physics->instance()->f;
			}
			if (_configuration.convergenceParameters().checkResidual()) {
				heatResidual = physics->sumSquares(physics->instance()->f, Physics::SumOperation::SUM, Physics::SumRestriction::NON_DIRICHLET, step.step);
				heatResidual += physics->sumSquares(physics->instance()->R, Physics::SumOperation::SUM, Physics::SumRestriction::DIRICHLET, step.step);
				heatResidual = sqrt(heatResidual);
				if (heatResidual < 1e-6) {
					heatResidual = 1e-6;
				}
			}
			sum(step, Matrices::f, Matrices::R, 1, -1);
			if (_configuration.convergenceParameters().checkResidual()) {
				sum(f_R_BtLambda, instance->f, instance->dualSolution, 1 , -1, "heat residual", "f - R", "Bt * Lambda");
				heatResidual = sqrt(physics->sumSquares(f_R_BtLambda, Physics::SumOperation::SUM)) / heatResidual;
				if (heatResidual < _configuration.convergenceParameters().requestedResidual()) {
					if (_configuration.convergenceParameters().checkSolution()) {
						if (temperatureResidual < _configuration.convergenceParameters().requestedSolution()) {
							break;
						}
					}
					break;
				}
			}

			composeGluing(step, Matrices::B1);
			multiply(step, Matrices::B1c, (step.iteration + 1) / (double)substeps);
			sum(step, Matrices::B1c, Matrices::primal, 1, -1);
			regularizeMatrices(step, Matrices::K);

			if (_configuration.method == NonLinearSolverBase::METHOD::MODIFIED_NEWTON_RHAPSON && step.substep) {
				updateLinearSolver(Matrices::f | Matrices::B1c);
			} else {
				updateLinearSolver(Matrices::K | Matrices::f | Matrices::B1c);
			}
			runLinearSolver();

			if (_configuration.line_search) {
				lineSearch(T, physics->instance()->primalSolution, F_ext, physics, step);
			}
			if (_configuration.convergenceParameters().checkSolution()) {
				temperatureResidual = sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));
			}
			sum(step, Matrices::primal, T, 1, 1, "prevSolution");
			if (_configuration.convergenceParameters().checkSolution()) {
				double temp = sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));
				if (temp < 1e-3) {
					temp = 1e-3;
				}
				temperatureResidual /= temp;
			}

			processSolution(step);
			storeSubSolution(step);
		}
	}
	step.iteration = substeps - 1;
}

void NewtonRhapson::postprocess(Step &step)
{
	processSolution(step);
	storeSolution(step);
}

void NewtonRhapson::finalize(Step &step)
{
	finalizeLinearSolver();
}



