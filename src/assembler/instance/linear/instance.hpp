
#include "instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void LinearInstance<TConstrains, TPhysics>::init()
{
	TimeEvent timePreparation("Prepare mesh structures"); timePreparation.start();
	_physics.prepareMeshStructures();
	timePreparation.endWithBarrier(); _timeStatistics.addEvent(timePreparation);


	TimeEvent timePhysics("Assemble stiffness matrices"); timePhysics.start();
	_physics.assembleStiffnessMatrices();
	timePhysics.endWithBarrier(); _timeStatistics.addEvent(timePhysics);

	if (config::info::PRINT_MATRICES) {
		_physics.save();
	}


	TimeEvent timeConstrains("Assemble gluing matrices"); timeConstrains.startWithBarrier();
	_physics.assembleGluingMatrices();
	timeConstrains.end(); _timeStatistics.addEvent(timeConstrains);

	if (config::info::PRINT_MATRICES) {
		_constrains.save();
	}


	TimeEvent timeSolver("Initialize solver"); timeSolver.startWithBarrier();
	_linearSolver.init(_mesh.neighbours());
	timeSolver.end(); _timeStatistics.addEvent(timeSolver);
}

template <class TConstrains, class TPhysics>
void LinearInstance<TConstrains, TPhysics>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	_linearSolver.Solve(_physics.f, solution);
	timeSolve.endWithBarrier(); _timeStatistics.addEvent(timeSolve);
}

template <class TConstrains, class TPhysics>
void LinearInstance<TConstrains, TPhysics>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
