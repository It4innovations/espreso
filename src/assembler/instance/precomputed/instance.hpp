
#include "instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void PrecomputedInstance<TConstrains, TPhysics>::init()
{
	TimeEvent timePreparation("Prepare mesh structures"); timePreparation.start();
	_physics.prepareMeshStructures();
	timePreparation.endWithBarrier(); _timeStatistics.addEvent(timePreparation);


	TimeEvent timePhysics("Assemble stiffness matrices"); timePhysics.start();
	_physics.assembleStiffnessMatrices();
	timePhysics.endWithBarrier(); _timeStatistics.addEvent(timePhysics);

	TimeEvent timeScaling("Assemble scaling matrices"); timeScaling.start();
	_physics.assembleScalingMatrices();
	timeScaling.endWithBarrier(); _timeStatistics.addEvent(timeScaling);

	if (config::info::PRINT_MATRICES) {
		_physics.saveStiffnessMatrices();
	}

	TimeEvent timeReg("Make K regular"); timeReg.start();
	_physics.makeStiffnessMatricesRegular();
	timeReg.endWithBarrier(); _timeStatistics.addEvent(timeReg);

	if (config::info::PRINT_MATRICES) {
		_physics.saveKernelMatrices();
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
void PrecomputedInstance<TConstrains, TPhysics>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();

	std::vector<std::vector<double> > tmpSolution(_physics.f.size());
	for (size_t i = 0; i <_physics.f.size(); i++) {
		tmpSolution[i].resize(_physics.f[i].size());
	}

	_linearSolver.Solve(_physics.f, tmpSolution);

	std::for_each(solution[0].begin(), solution[0].end(), [] (double &v) { v = 0; });

	for (size_t p = 0; p < _mesh.parts(); p++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(p);
		for (size_t i = 0; i < l2c.size(); i++) {
			for (size_t d = 0; d < _physics.pointDOFs.size(); d++) {
				solution[0][_physics.pointDOFs.size() * l2c[i] + d] += tmpSolution[p][_physics.pointDOFs.size() * i + d] / _mesh.nodes()[l2c[i]]->domains().size();
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
