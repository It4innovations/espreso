
#include "instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void LinearInstance<TConstrains, TPhysics>::init()
{
	TimeEvent timePreparation("Prepare mesh structures"); timePreparation.start();
	_physics.prepareMeshStructures();
	timePreparation.endWithBarrier(); _timeStatistics.addEvent(timePreparation);

	if (config::output::SAVE_PROPERTIES || config::output::SAVE_RESULTS) {
		_store.storeGeometry();
	}
	if (config::output::SAVE_PROPERTIES) {
		_physics.saveMeshProperties(_store);
	}

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

	if (config::output::SAVE_GLUING) {
		output::VTK::gluing(_mesh, _constrains, "B1", _physics.pointDOFs.size(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
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

	if (config::output::SAVE_RESULTS) {
		_physics.saveMeshResults(_store, solution);
	}
}

template <class TConstrains, class TPhysics>
void LinearInstance<TConstrains, TPhysics>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
