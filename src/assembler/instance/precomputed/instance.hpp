
#include "instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void PrecomputedInstance<TConstrains, TPhysics>::init()
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
	_linearSolver.init(_mesh.neighbours());
	timeSolver.end(); _timeStatistics.addEvent(timeSolver);
}

template <class TConstrains, class TPhysics>
void PrecomputedInstance<TConstrains, TPhysics>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();

	std::vector<std::vector<double> > tmpSolution(_physics.f.size());
	for (size_t i = 0; i <_physics.f.size(); i++) {
		tmpSolution[i].resize(_physics.f[i].size());
	}

	_linearSolver.Solve(_physics.f, tmpSolution);

	std::for_each(solution[0].begin(), solution[0].end(), [] (double &v) { v = 0; });

	const Boundaries &boundaries = _mesh.subdomainBoundaries();
	for (size_t p = 0; p < _mesh.parts(); p++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(p);
		for (size_t i = 0; i < l2c.size(); i++) {
			for (size_t d = 0; d < _physics.DOFs.size(); d++) {
				solution[0][_physics.DOFs.size() * l2c[i] + d] += tmpSolution[p][_physics.DOFs.size() * i + d] / boundaries[l2c[i]].size();
			}
		}
	}

	timeSolve.endWithBarrier(); _timeStatistics.addEvent(timeSolve);
}

template <class TConstrains, class TPhysics>
void PrecomputedInstance<TConstrains, TPhysics>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
