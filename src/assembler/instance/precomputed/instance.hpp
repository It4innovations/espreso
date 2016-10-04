
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

	const APIMesh &mesh = static_cast<const APIMesh&>(_mesh);

	size_t DOFIndex = 0;
	for (size_t i = 0; i < mesh.DOFs().size(); i++) {
		solution[0][i] = 0;
		for (size_t d = 0; d < mesh.DOFs()[i]->domains().size(); d++) {
			solution[0][i] += tmpSolution[mesh.DOFs()[i]->domains()[d]][mesh.DOFs()[i]->DOFIndex(mesh.DOFs()[i]->domains()[d], DOFIndex)];
		}
		solution[0][i] /= mesh.DOFs()[i]->domains().size();
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
