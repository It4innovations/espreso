
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
		output::Store* store,
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
	multiply(step, Matrices::f, (step.substep + 1) / (double)substeps);
	composeGluing(step, Matrices::B1);
	multiply(step, Matrices::B1c, (step.substep + 1) / (double)substeps);
	regularizeMatrices(step, Matrices::K);
	composeGluing(step, Matrices::B0);
}

void NewtonRhapson::preprocess(Step &step)
{
	initLinearSolver(step);
}

void NewtonRhapson::solve(Step &step)
{
	if (!_configuration.convergenceParameters().checkSolution() && !_configuration.convergenceParameters().checkResidual()) {
		ESINFO(GLOBAL_ERROR) << "It is not possible to turn off the both 'temperature' and 'heat' convergence.";
	}

	runLinearSolver(step);
	processSolution(step);
	storeSubSolution(step);

	std::vector<std::vector<double> > T;
	std::vector<std::vector<double> > F_ext;
	std::vector<std::vector<double> > f_R_BtLambda;

	double alpha, max;
	int cumiter = 0;

	size_t substeps = _configuration.stepping == NonLinearSolverBase::STEPPINGG::TRUE ? _configuration.substeps : 1;
	for (step.substep = 0; step.substep < substeps; step.substep++) {

		ESINFO(CONVERGENCE) <<  " ";
		ESINFO(CONVERGENCE) <<  "**** LOAD_STEP " << step.step + 1 << " SUBSTEP "<< step.substep + 1 << " START";

		step.iteration = 0;
		double temperatureResidual = 10 * _configuration.convergenceParameters().requestedSolution();
		double temperatureResidual_first = 0;
		double temperatureResidual_second = 0;

		double heatResidual;
		double heatResidual_first = 0;
		double heatResidual_second = 0;

		while (
			step.iteration++ < _configuration.max_iterations ){

			cumiter +=1;
			if (!_configuration.convergenceParameters().checkResidual()) {
				ESINFO(CONVERGENCE) <<  " ";
				ESINFO(CONVERGENCE) <<  " >> EQUILIBRIUM ITERATION " << step.iteration + 1 << " IN SUBSTEP "  << step.substep + 1;
			}

			T = physics->instance()->primalSolution;

			if (_configuration.method == NonLinearSolverBase::METHOD::MODIFIED_NEWTON_RHAPSON && step.iteration > 1) {
				updateMatrices(step, Matrices::f | Matrices::R, physics->instance()->solutions);
			} else {
				updateMatrices(step, Matrices::K | Matrices::f | Matrices::R, physics->instance()->solutions);
			}
			multiply(step, Matrices::f, (step.substep + 1) / (double)substeps);

			if (_configuration.line_search) {
				F_ext = physics->instance()->f;
			}
			if (_configuration.convergenceParameters().checkResidual()) {
				heatResidual_second = physics->sumSquares(physics->instance()->f, Physics::SumOperation::SUM, Physics::SumRestriction::NON_DIRICHLET, step.step);
				heatResidual_second += physics->sumSquares(physics->instance()->R, Physics::SumOperation::SUM, Physics::SumRestriction::DIRICHLET, step.step);
				heatResidual_second = sqrt(heatResidual_second);
				if (heatResidual_second < 1e-3) {
					heatResidual_second = 1e-3;
				}
			}
			sum(step, Matrices::f, Matrices::R, 1, -1);
			if (_configuration.convergenceParameters().checkResidual()) {
				sum(f_R_BtLambda, instance->f, instance->dualSolution, 1 , -1, "heat residual", "f - R", "Bt * Lambda");
				heatResidual_first = sqrt(physics->sumSquares(f_R_BtLambda, Physics::SumOperation::SUM));
				heatResidual = heatResidual_first / heatResidual_second;

				if (heatResidual < _configuration.convergenceParameters().requestedResidual()  && step.iteration > 1 ) {

					ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  heatResidual_first << "  CRITERION_VALUE = " << heatResidual_second * _configuration.convergenceParameters().requestedResidual() << " <<< CONVERGED >>>";

					if (_configuration.convergenceParameters().checkSolution()) {
						if (temperatureResidual < _configuration.convergenceParameters().requestedSolution()) {
							break;
						}
					}else{
						break;
					}

				}else{

					ESINFO(CONVERGENCE) <<  " ";
					ESINFO(CONVERGENCE) <<  " >> EQUILIBRIUM ITERATION " << step.iteration + 1 << " IN SUBSTEP "  << step.substep + 1;
					ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  heatResidual_first << "  CRITERION_VALUE = " << heatResidual_second * _configuration.convergenceParameters().requestedResidual();
				}
			}

			composeGluing(step, Matrices::B1);
			multiply(step, Matrices::B1c, (step.substep + 1) / (double)substeps);
			sum(step, Matrices::B1c, Matrices::primal, 1, -1);
			regularizeMatrices(step, Matrices::K);

			if (_configuration.method == NonLinearSolverBase::METHOD::MODIFIED_NEWTON_RHAPSON && step.iteration) {
				updateLinearSolver(step, Matrices::f | Matrices::B1c);
			} else {
				updateLinearSolver(step, Matrices::K | Matrices::f | Matrices::B1c);
			}
			runLinearSolver(step);

			ESINFO(CONVERGENCE) <<  "    LINEAR_SOLVER_OUTPUT: SOLVER = " << "PCG" <<   " N_ITERATIONS = " << "1" << "  " ;



			if (_configuration.line_search && step.iteration > 1) {
				max = maxAbsValue(physics->instance()->primalSolution);
				alpha = lineSearch(T, physics->instance()->primalSolution, F_ext, physics, step);
				ESINFO(CONVERGENCE) << "    LINE_SEARCH_OUTPUT: " << "PARAMETER = " << alpha << "  MAX_DOF_INCREMENT = " << max << "  SCALED_MAX_INCREMENT = "  <<  alpha*max ;
			} else if (_configuration.line_search && step.iteration == 1){
				max = maxAbsValue(physics->instance()->primalSolution);
				ESINFO(CONVERGENCE) << "    LINE_SEARCH_OUTPUT: " << "PARAMETER = 1" << "  MAX_DOF_INCREMENT = " << max << "  SCALED_MAX_INCREMENT = "  <<  max ;
			}
			if (_configuration.convergenceParameters().checkSolution()) {
				temperatureResidual_first = sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));
			}
			sum(step, Matrices::primal, T, 1, 1, "prevSolution");
			if (_configuration.convergenceParameters().checkSolution()) {
				 temperatureResidual_second = sqrt(physics->sumSquares(physics->instance()->primalSolution, Physics::SumOperation::AVERAGE));
				if (temperatureResidual_second < 1e-3) {
					temperatureResidual_second = 1e-3;
				}
				temperatureResidual = temperatureResidual_first / temperatureResidual_second;

				if ( temperatureResidual > _configuration.convergenceParameters().requestedSolution()){
					ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  temperatureResidual_first << "  CRITERION_VALUE = " << temperatureResidual_second * _configuration.convergenceParameters().requestedSolution() ;
				} else {
					ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  temperatureResidual_first << "  CRITERION_VALUE = " << temperatureResidual_second * _configuration.convergenceParameters().requestedSolution() <<  " <<< CONVERGED >>>" ;
					if (!_configuration.convergenceParameters().checkResidual()){
						break;
					}
				}
			}


			processSolution(step);
			storeSubSolution(step);
		}
		if (_configuration.convergenceParameters().checkResidual()) {
			ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << step.iteration ;
		}else{
			ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << step.iteration + 1 ;
		}
		ESINFO(CONVERGENCE) <<  " >> SUBSTEP " << step.substep + 1 << " IS DONE /" << " CUMULATIVE ITERATION NUMBER = " << cumiter + 1 ;
		ESINFO(CONVERGENCE) <<  " ";

	}
	step.substep = substeps - 1;
}

void NewtonRhapson::postprocess(Step &step)
{
	processSolution(step);
	storeSolution(step);
}

void NewtonRhapson::finalize(Step &step)
{
	finalizeLinearSolver(step);
}



