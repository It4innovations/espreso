
#include "loadstepsolver.h"
#include "../timestep/timestepsolver.h"

#include "../../step.h"
#include "../assembler.h"

#include "../../../basis/logging/logging.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(const std::string &description, TimeStepSolver &timeStepSolver, double duration)
: _description(description), _timeStepSolver(timeStepSolver), _assembler(timeStepSolver._assembler), _duration(duration), _startTime(0), _precision(1e-6)
{

}

std::string LoadStepSolver::description() const
{
	return _description;
}

void LoadStepSolver::initLoadStep(Step &step)
{
	if (step.step == 0) {
		_assembler.preprocessData(step);
	}
	_assembler.setRegularizationCallback();
	_assembler.setB0Callback();
}

bool LoadStepSolver::hasNextTimeStep(Step &step)
{
	return step.currentTime + _precision < _startTime + _duration;
}

void LoadStepSolver::finalizeLoadStep(Step &step)
{
	_assembler.finalize();
}

void LoadStepSolver::run(Step &step)
{
	ESINFO(PROGRESS1) << "Solve LOAD STEP " << step.step + 1 << ": " << description() << " with " << _timeStepSolver.description() << " time step(s).";

	_startTime = step.currentTime;
	step.substep = 0;
	step.iteration = 0;

	initLoadStep(step);
	while (hasNextTimeStep(step)) {
		runNextTimeStep(step);
		ESINFO(PROGRESS1) << description() << " SOLVER: load step " << step.step << ", time step " << step.substep + 1 << "[" << step.currentTime << "s] finished.";
		step.substep++;
		step.iteration = 0;
	}
	finalizeLoadStep(step);
}
