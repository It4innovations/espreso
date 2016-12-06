
#include "instance.h"

namespace espreso {

template <class TPhysics>
void SemiSmoothNewtonMethod<TPhysics>::init()
{
	TimeEvent timePreparation("Prepare mesh structures"); timePreparation.start();
	_physics.prepareMeshStructures();
	timePreparation.endWithBarrier(); _timeStatistics.addEvent(timePreparation);

	if (output->properties || output->results) {
		_store.storeGeometry();
	}
	if (output->properties) {
		_physics.saveMeshProperties(_store);
	}

	TimeEvent timePhysics("Assemble stiffness matrices"); timePhysics.start();
	_physics.assembleStiffnessMatrices();
	timePhysics.endWithBarrier(); _timeStatistics.addEvent(timePhysics);

	if (output->print_matrices) {
		_physics.saveStiffnessMatrices();
	}

	TimeEvent timeReg("Make K regular"); timeReg.start();
	_physics.makeStiffnessMatricesRegular();
	timeReg.endWithBarrier(); _timeStatistics.addEvent(timeReg);

	if (output->print_matrices) {
		_physics.saveKernelMatrices();
	}

	TimeEvent timeConstrains("Assemble gluing matrices"); timeConstrains.startWithBarrier();
	_physics.assembleB1();
	timeConstrains.end(); _timeStatistics.addEvent(timeConstrains);

	if (output->print_matrices) {
		_constrains.save();
	}

	if (output->gluing) {
		store::VTK::gluing(_mesh, _constrains, "B1", _physics.pointDOFs.size(), output->domain_shrink_ratio, output->cluster_shrink_ratio);
	}

	TimeEvent timeSolver("Initialize solver"); timeSolver.startWithBarrier();
	_linearSolver.init(_mesh.neighbours());
	timeSolver.end(); _timeStatistics.addEvent(timeSolver);
}

template <class TPhysics>
void SemiSmoothNewtonMethod<TPhysics>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	_linearSolver.Solve(_physics.f, solution);
	timeSolve.endWithBarrier(); _timeStatistics.addEvent(timeSolve);

	if (output->results) {
		_physics.saveMeshResults(_store, solution);
	}
}

template <class TPhysics>
void SemiSmoothNewtonMethod<TPhysics>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
