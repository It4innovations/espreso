
#include "instance.h"

namespace espreso {

template <class TPhysics, class TConfiguration>
void LinearInstance<TPhysics, TConfiguration>::init()
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

	TimeEvent timeConstrainsB1("Assemble B1"); timeConstrainsB1.startWithBarrier();
	_physics.assembleB1();
	timeConstrainsB1.end(); _timeStatistics.addEvent(timeConstrainsB1);

	TimeEvent timeReg("Make K regular"); timeReg.start();
	_physics.makeStiffnessMatricesRegular();
	timeReg.endWithBarrier(); _timeStatistics.addEvent(timeReg);

	if (output->print_matrices) {
		_physics.saveKernelMatrices();
	}

	TimeEvent timeConstrainsB0("Assemble B0"); timeConstrainsB0.startWithBarrier();
	_physics.assembleB0();
	timeConstrainsB0.end(); _timeStatistics.addEvent(timeConstrainsB0);

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

template <class TPhysics, class TConfiguration>
void LinearInstance<TPhysics, TConfiguration>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	_linearSolver.Solve(_physics.f, solution);
	timeSolve.endWithBarrier(); _timeStatistics.addEvent(timeSolve);

	if (output->results) {
		_physics.saveMeshResults(_store, solution);
	}
}

template <class TPhysics, class TConfiguration>
void LinearInstance<TPhysics, TConfiguration>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
