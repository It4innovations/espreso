
#include "pseudotimestepping.h"

#include "../assembler.h"
#include "../timestep/timestepsolver.h"
#include "../../step.h"
#include "../../instance.h"
#include "../../../configuration/physics/nonlinearsolver.h"


using namespace espreso;

PseudoTimeStepping::PseudoTimeStepping(TimeStepSolver &timeStepSolver, const NonLinearSolverBase &configuration, double duration)
: LoadStepSolver("PSEUDO TIME STEPS", timeStepSolver, duration), _configuration(configuration)
{
	_assembler.setRegularizationCallback();
	_assembler.setB0Callback();
}

Matrices PseudoTimeStepping::updateStructuralMatrices(Step &step, Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::f | Matrices::R | Matrices::B1);

//	TODO: ??
//	if (step.substep && !_timeDependent) {
//		updatedMatrices &= (Matrices::f | Matrices::B1); // TODO: update only B1c
//	}

	return reassembleStructuralMatrices(step, updatedMatrices);
}

Matrices PseudoTimeStepping::reassembleStructuralMatrices(Step &step, Matrices matrices)
{
	_assembler.updateMatrices(step, matrices);
	return matrices;
}

void PseudoTimeStepping::runNextTimeStep(Step &step)
{
	double last = step.currentTime;
	step.currentTime += _duration / _configuration.substeps;
	if (step.currentTime + _precision >= _startTime + _duration) {
		step.currentTime = _startTime + _duration;
	}
	step.timeStep = step.currentTime - last;
	processTimeStep(step);
}

void PseudoTimeStepping::processTimeStep(Step &step)
{
	step.internalForceReduction = (double)(step.substep + 1) / _configuration.substeps;
	step.timeIntegrationConstantK = 1;
	step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(step, *this);

	_assembler.processSolution(step);
	_assembler.storeSolution(step);
}




