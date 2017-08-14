
#include "newtonrhapson.h"

#include "../assembler.h"
#include "../loadstep/loadstepsolver.h"
#include "../../physics/physics.h"
#include "../../instance.h"
#include "../../step.h"

#include "../../../configuration/physics/nonlinearsolver.h"
#include "../../../basis/logging/logging.h"
#include "../../../linearsolver/linearsolver.h"

#include <cmath>

using namespace espreso;

NewtonRhapson::NewtonRhapson(Assembler &assembler, const NonLinearSolverBase &configuration)
: TimeStepSolver("Newton Rhapson", assembler), _configuration(configuration)
{

}

void NewtonRhapson::solve(Step &step, LoadStepSolver &loadStepSolver)
{
	if (!_configuration.convergenceParameters().checkSolution() && !_configuration.convergenceParameters().checkResidual()) {
		ESINFO(GLOBAL_ERROR) << "Turn on at least one convergence parameter for NONLINEAR solver.";
	}

	Matrices updatedMatrices;
	double &solverPrecision = _assembler.linearSolver.precision();
	double solverPrecisionError = 1;

	double temperatureResidual = 10 * _configuration.convergenceParameters().requestedSolution();
	double temperatureResidual_first = 0;
	double temperatureResidual_second = 0;

	double heatResidual;
	double heatResidual_first = 0;
	double heatResidual_second = 0;

	double alpha, maxSolutionValue;


	step.iteration = 0;
	_assembler.solve(step, loadStepSolver.updateStructuralMatrices(step, Matrices::K | Matrices::M | Matrices::f | Matrices::B1));
	_assembler.processSolution(step);
	_assembler.storeSubSolution(step);

	while (step.iteration++ < _configuration.max_iterations) {
		if (!_configuration.convergenceParameters().checkResidual()) {
			ESINFO(CONVERGENCE) << "\n >> EQUILIBRIUM ITERATION " << step.iteration + 1 << " IN SUBSTEP "  << step.substep + 1;
		}

		_solution = _assembler.instance.primalSolution;
		if (_configuration.method == NonLinearSolverBase::METHOD::NEWTON_RHAPSON) {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(step, Matrices::K | Matrices::M | Matrices::f | Matrices::R);
		} else {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(step, Matrices::f | Matrices::R);
		}
		if (_configuration.line_search) {
			_f_ext = _assembler.instance.f;
		}
		if (_configuration.convergenceParameters().checkResidual()) {
			heatResidual_second = _assembler.sumSquares(step, _assembler.instance.f, SumOperation::SUM, SumRestriction::NON_DIRICHLET, "norm of f not on DIRICHLET");
			heatResidual_second += _assembler.sumSquares(step, _assembler.instance.R, SumOperation::SUM, SumRestriction::DIRICHLET, "norm of R on DIRICHLET");
			heatResidual_second = sqrt(heatResidual_second);
			if (heatResidual_second < 1e-3) {
				heatResidual_second = 1e-3;
			}
		}

		_assembler.sum(
				_assembler.instance.f,
				1, _assembler.instance.f,
				-1, _assembler.instance.R,
				"f = f - R");

		if (_configuration.convergenceParameters().checkResidual()) {
			_assembler.sum(
					_f_R_BtLambda,
					1, _assembler.instance.f,
					-1, _assembler.instance.dualSolution,
					"(f - R) * Bt * Lambda");

			heatResidual_first = sqrt(_assembler.sumSquares(step, _f_R_BtLambda, SumOperation::SUM, SumRestriction::NONE, "norm of (f - R) * Bt * Lambda"));
			heatResidual = heatResidual_first / heatResidual_second;

			if (heatResidual < _configuration.convergenceParameters().requestedResidual() && step.iteration > 1 ) {
				ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  heatResidual_first << "  CRITERION_VALUE = " << heatResidual_second * _configuration.convergenceParameters().requestedResidual() << " <<< CONVERGED >>>";
				if (_configuration.convergenceParameters().checkSolution()) {
					if (temperatureResidual < _configuration.convergenceParameters().requestedSolution()) {
						break;
					}
				} else {
					break;
				}
			} else {
				ESINFO(CONVERGENCE) <<  "]n >> EQUILIBRIUM ITERATION " << step.iteration + 1 << " IN SUBSTEP "  << step.substep + 1;
				ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  heatResidual_first << "  CRITERION_VALUE = " << heatResidual_second * _configuration.convergenceParameters().requestedResidual();
			}
		}

		if (updatedMatrices & Matrices::K) {
			updatedMatrices |= loadStepSolver.reassembleStructuralMatrices(step, Matrices::B1c | Matrices::B1duplicity);
		} else {
			updatedMatrices |= loadStepSolver.reassembleStructuralMatrices(step, Matrices::B1c);
		}
		_assembler.addToDirichletInB1(-1, _assembler.instance.primalSolution);

		if (_configuration.adaptive_precision) {
			if (step.iteration > 1) {
				solverPrecisionError = temperatureResidual_first / temperatureResidual_second;
				solverPrecision = std::min(_configuration.r_tol * solverPrecisionError, _configuration.c_fact * solverPrecision);
			}
			ESINFO(CONVERGENCE) << "    ADAPTIVE PRECISION = " << solverPrecision << " EPS_ERR = " << solverPrecisionError;
		}

		_assembler.solve(step, updatedMatrices);
		ESINFO(CONVERGENCE) <<  "    LINEAR_SOLVER_OUTPUT: SOLVER = " << "PCG" <<   " N_ITERATIONS = " << "1" << "  " ;

		if (_configuration.line_search) {
			maxSolutionValue =_assembler.maxAbsValue(_assembler.instance.primalSolution, "max = |solution|");
			alpha = _assembler.lineSearch(step, _solution, _assembler.instance.primalSolution, _f_ext);
			ESINFO(CONVERGENCE) << "    LINE_SEARCH_OUTPUT: " << "PARAMETER = " << alpha << "  MAX_DOF_INCREMENT = " << maxSolutionValue << "  SCALED_MAX_INCREMENT = " << alpha * maxSolutionValue;
		}
		if (_configuration.convergenceParameters().checkSolution()) {
			temperatureResidual_first = sqrt(_assembler.sumSquares(step, _assembler.instance.primalSolution, SumOperation::AVERAGE, SumRestriction::NONE, "|delta U|"));
		}
		_assembler.sum(
				_assembler.instance.primalSolution,
				1, _assembler.instance.primalSolution,
				1, _solution, "U = delta U + U");

		if (_configuration.convergenceParameters().checkSolution()) {
			 temperatureResidual_second = sqrt(_assembler.sumSquares(step, _assembler.instance.primalSolution, SumOperation::AVERAGE, SumRestriction::NONE, "|U|"));
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

		_assembler.processSolution(step);
		_assembler.storeSubSolution(step);
	}

	if (_configuration.convergenceParameters().checkResidual()) {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << step.iteration ;
	} else {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << step.iteration + 1 ;
	}
}


