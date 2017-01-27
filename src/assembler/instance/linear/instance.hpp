
#include "instance.h"

#include "../../../mesh/structures/mesh.h"

namespace espreso {

template <class TPhysics, class TConfiguration>
void LinearInstance<TPhysics, TConfiguration>::init()
{
	TimeEvent timePreparation("Prepare mesh structures"); timePreparation.start();
	_physics.prepareMeshStructures();
	timePreparation.endWithBarrier(); _timeStatistics.addEvent(timePreparation);

	if (_output.properties || _output.results) {
		_store.storeGeometry();
	}
	if (_output.properties) {
		_physics.saveMeshProperties(_store);
	}

	TimeEvent timePhysics("Assemble stiffness matrices"); timePhysics.start();
	_physics.assembleStiffnessMatrices();
	timePhysics.endWithBarrier(); _timeStatistics.addEvent(timePhysics);

	if (environment->print_matrices) {
		_physics.saveStiffnessMatrices();
	}

	TimeEvent timeConstrainsB1("Assemble B1"); timeConstrainsB1.startWithBarrier();
	_physics.assembleB1();
	timeConstrainsB1.end(); _timeStatistics.addEvent(timeConstrainsB1);

	TimeEvent timeReg("Make K regular"); timeReg.start();
	_physics.makeStiffnessMatricesRegular();
	timeReg.endWithBarrier(); _timeStatistics.addEvent(timeReg);

	if (environment->print_matrices) {
		_physics.saveKernelMatrices();
	}

	TimeEvent timeConstrainsB0("Assemble B0"); timeConstrainsB0.startWithBarrier();
	_physics.assembleB0();
	timeConstrainsB0.end(); _timeStatistics.addEvent(timeConstrainsB0);

	if (environment->print_matrices) {
		_constrains.save();
	}

	if (_output.gluing) {
		store::VTK::gluing(_output, _mesh, _constrains, "B1", _physics.pointDOFs.size());
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
	if (_output.results) {
		_physics.postProcess(_store, solution);
	}
	timeSolve.endWithBarrier(); _timeStatistics.addEvent(timeSolve);

	if (_output.results) {
		_physics.saveMeshResults(_store, solution);
	}
	if (_output.properties || _output.results) {
		_store.finalize();
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
