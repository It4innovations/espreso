
#include "steadystate.h"
#include "../timestep/timestepsolver.h"

#include "../../step.h"
#include "../../instance.h"
#include "../assembler.h"

using namespace espreso;

SteadyStateSolver::SteadyStateSolver(TimeStepSolver &timeStepSolver, double duration)
: LoadStepSolver("STEADY STATE", timeStepSolver, duration)
{
	_assembler.setRegularizationCallback();
	_assembler.setB0Callback();
}

Matrices SteadyStateSolver::updateStructuralMatrices(Step &step, Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::f | Matrices::R | Matrices::B1 | Matrices::B1c | Matrices::B1duplicity);

	return reassembleStructuralMatrices(step, updatedMatrices);
}

Matrices SteadyStateSolver::reassembleStructuralMatrices(Step &step, Matrices matrices)
{
	_assembler.updateMatrices(step, matrices);
	return matrices;
}

void SteadyStateSolver::runNextTimeStep(Step &step)
{
	step.currentTime += _duration;
	step.timeStep = _duration;
	processTimeStep(step);
}


void SteadyStateSolver::processTimeStep(Step &step)
{
	step.internalForceReduction = 1;
	step.timeIntegrationConstantK = 1;
	step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(step, *this);

	_assembler.processSolution(step);
	_assembler.storeSolution(step);
}

