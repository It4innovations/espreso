
#include "instance.h"

namespace espreso {

template <class TPhysics>
void DynamicsInstance<TPhysics>::init()
{
	TimeEvent timePreparation("Prepare mesh structures"); timePreparation.start();
	_physics.prepareMeshStructures();
	timePreparation.endWithBarrier(); _timeStatistics.addEvent(timePreparation);

	if (output->properties || output->results) {
		_store.storeGeometry(_time);
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


	TimeEvent timeSolver("Initialize solver"); timeSolver.startWithBarrier();
	_linearSolver.init(_mesh.neighbours());
	timeSolver.end(); _timeStatistics.addEvent(timeSolver);

	TimeEvent timePrep("Prepare vectors for transient problem"); timePrep.startWithBarrier();

	 _u.resize(_mesh.parts());
	 _v.resize(_mesh.parts());
	 _w.resize(_mesh.parts());

	 _vn.resize(_mesh.parts());
	 _wn.resize(_mesh.parts());

	 _b.resize(_mesh.parts());
	 _tmp.resize(_mesh.parts());

	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		_u[p].resize(_physics.M[p].rows, 0);
		_v[p].resize(_physics.M[p].rows, 0);
		_w[p].resize(_physics.M[p].rows, 0);

		_vn[p].resize(_physics.M[p].rows, 0);
		_wn[p].resize(_physics.M[p].rows, 0);

		_b[p].resize(_physics.M[p].rows, 0);
		_tmp[p].resize(_physics.M[p].rows, 0);
	}

	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t i = 1; i < _w[p].size(); i += 3) {
			_w[p][i] = 1.0;
		}
	}

	timePrep.end(); _timeStatistics.addEvent(timePrep);
}

template <class TPhysics>
void DynamicsInstance<TPhysics>::solve(std::vector<std::vector<double> > &solution)
{
	ESINFO(PROGRESS1) << "Time: " << _time;

	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		for(size_t i = 0; i < _u[p].size(); i++) {
			_tmp[p][i] = _physics.A[0] * _u[p][i] + _physics.A[2] * _v[p][i] + _physics.A[3] * _w[p][i];
		}
		_physics.M[p].MatVec(_tmp[p], _b[p], 'N');
	}

	if (_time && output->results) {
		_store.storeGeometry(_time);
	}

	TimeEvent timeLSrun("Linear Solver - runtime"); timeLSrun.start();
	solution.resize(_mesh.parts());
	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		solution[p].resize(_physics.M[p].rows, 0);
	}

	_linearSolver.Solve(_b, solution);
	timeLSrun.endWithBarrier(); _timeStatistics.addEvent(timeLSrun);

	if (output->results) {
		_physics.saveMeshResults(_store, solution);
	}

	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		for(size_t i = 0; i < _u[p].size(); i++) {
			_wn[p][i] = (_physics.A[0] * (solution[p][i] - _u[p][i])) - (_physics.A[2] * _v[p][i]) - (_physics.A[3] * _w[p][i]);
			_vn[p][i] = _v[p][i] + (_physics.A[6] * _w[p][i]) + (_physics.A[7] * _wn[p][i]);

			_u[p][i] = solution[p][i];
			_v[p][i] = _vn[p][i];
			_w[p][i] = _wn[p][i];
		}
	}

	_time++;
}

template <class TPhysics>
void DynamicsInstance<TPhysics>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
