
#include "instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void LinearInstance<TConstrains, TPhysics>::init()
{
	TimeEvent timePhysics("Assemble problem"); timePhysics.start();
	ESINFO(PROGRESS2) << "Assemble problem";
	_physics.assemble();
	timePhysics.endWithBarrier(); _timeStatistics.addEvent(timePhysics);
	if (config::info::PRINT_MATRICES) {
		_physics.save();
	}


	TimeEvent timeConstrains("Equality constrains"); timeConstrains.startWithBarrier();
	ESINFO(PROGRESS2) << "Assemble equality constraints";
	_constrains.assemble();
	timeConstrains.end(); _timeStatistics.addEvent(timeConstrains);
	if (config::info::PRINT_MATRICES) {
		_constrains.save();
	}


	TimeEvent timeSolver("Initialize solver"); timeSolver.startWithBarrier();
	ESINFO(PROGRESS2) << "Initialize linear solver";
	_linearSolver.DOFS_PER_NODE = _physics.DOFs;
	_linearSolver.setup(config::env::MPIrank, config::env::MPIsize, true);
	_linearSolver.init(
			_mesh,
			_physics.K,
			_physics.T,
			_constrains.B1,
			_constrains.B0,
			_constrains.B1subdomainsMap,
			_constrains.B0subdomainsMap,
			_constrains.B1clustersMap,
			_constrains.B1duplicity,
			_physics.f,
			_constrains.B1c,
			_mesh.getFixPoints(),
			_mesh.neighbours());
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
