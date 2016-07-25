
#include "instance.h"
#include "esoutput.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void DynamicsInstance<TConstrains, TPhysics>::init()
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
	_linearSolver.setup(config::env::MPIrank, config::env::MPIsize, false);
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


	TimeEvent timePrep("Prepare vectors for transient problem"); timePrep.startWithBarrier();

	 _u.resize(_mesh.parts());
	 _v.resize(_mesh.parts());
	 _w.resize(_mesh.parts());

	 _vn.resize(_mesh.parts());
	 _wn.resize(_mesh.parts());

	 _b.resize(_mesh.parts());
	 _tmp.resize(_mesh.parts());

	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		_u[p].resize(_physics.M[p].rows, 0);
		_v[p].resize(_physics.M[p].rows, 0);
		_w[p].resize(_physics.M[p].rows, 0);

		_vn[p].resize(_physics.M[p].rows, 0);
		_wn[p].resize(_physics.M[p].rows, 0);

		_b[p].resize(_physics.M[p].rows, 0);
		_tmp[p].resize(_physics.M[p].rows, 0);
	}

	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t i = 1; i < _w[p].size(); i += 3) {
			_w[p][i] = 1.0;
		}
	}

	timePrep.end(); _timeStatistics.addEvent(timePrep);
}

template <class TConstrains, class TPhysics>
void DynamicsInstance<TConstrains, TPhysics>::pre_solve_update(std::vector<std::vector<double> > &solution)
{
	ESINFO(PROGRESS1) << "Time: " << _time;

	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		for(size_t i = 0; i < _u[p].size(); i++) {
			_tmp[p][i] = _physics.A[0] * _u[p][i] + _physics.A[2] * _v[p][i] + _physics.A[3] * _w[p][i];
		}
		_physics.M[p].MatVec(_tmp[p], _b[p], 'N');
	}
}

template <class TConstrains, class TPhysics>
void DynamicsInstance<TConstrains, TPhysics>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeLSrun("Linear Solver - runtime"); timeLSrun.start();
	solution.resize(_mesh.parts());
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		solution[p].resize(_physics.M[p].rows, 0);
	}

	_linearSolver.Solve(_b, solution);
	timeLSrun.endWithBarrier(); _timeStatistics.addEvent(timeLSrun);
}

template <class TConstrains, class TPhysics>
void DynamicsInstance<TConstrains, TPhysics>::post_solve_update(std::vector<std::vector<double> > &solution)
{
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		for(size_t i = 0; i < _u[p].size(); i++) {
			_wn[p][i] = (_physics.A[0] * (solution[p][i] - _u[p][i])) - (_physics.A[2] * _v[p][i]) - (_physics.A[3] * _w[p][i]);
			_vn[p][i] = _v[p][i] + (_physics.A[6] * _w[p][i]) + (_physics.A[7] * _wn[p][i]);

			_u[p][i] = solution[p][i];
			_v[p][i] = _vn[p][i];
			_w[p][i] = _wn[p][i];
		}
	}

//#ifdef CATALYST
//	unsigned int timeStep = tt;
//	double time = timeStep * dynamic_timestep;
//	Adaptor::CoProcess(input.mesh,l2g_vec, vec_u_n,  time, timeStep, timeStep == numberOfTimeSteps - 1);
//#endif



//
//	_instance.mesh().store(mesh::VTK_FULL, ss.str(), vec_u_n, 0.95, 0.9);
//	saveVTK(ss.str().c_str(), vec_u_n, l2g_vec, _instance.localBoundaries(), _instance.globalBoundaries(), 0.95, 0.9);
//
//


	_time++;
}

template <class TConstrains, class TPhysics>
void DynamicsInstance<TConstrains, TPhysics>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
