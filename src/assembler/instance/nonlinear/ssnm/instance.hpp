
#include "instance.h"

namespace espreso {

template <class TPhysics, class TConfiguration>
void SemiSmoothNewtonMethod<TPhysics, TConfiguration>::init()
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

	TimeEvent timeReg("Make K regular"); timeReg.start();
	_physics.makeStiffnessMatricesRegular();
	timeReg.endWithBarrier(); _timeStatistics.addEvent(timeReg);

	if (environment->print_matrices) {
		_physics.saveKernelMatrices();
	}

	TimeEvent timeConstrains("Assemble gluing matrices"); timeConstrains.startWithBarrier();
	_physics.assembleB1();
	timeConstrains.end(); _timeStatistics.addEvent(timeConstrains);

	if (environment->print_matrices) {
		_constraints.save();
	}

	if (_output.gluing) {
		store::VTK::gluing(_output, _mesh, _constraints, "B1", _physics.pointDOFs.size());
	}

	TimeEvent timeSolver("Initialize solver"); timeSolver.startWithBarrier();
	_linearSolver.init(_mesh.neighbours());
	timeSolver.end(); _timeStatistics.addEvent(timeSolver);
}

template <class TPhysics, class TConfiguration>
void SemiSmoothNewtonMethod<TPhysics, TConfiguration>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeSolve("Linear Solver - runtime");

	auto computeError = [] (const Mesh &mesh, const std::vector<std::vector<double> > &prev, const std::vector<std::vector<double> > &current) {
		// norm(u - prev) / norm (u);
		std::vector<double> n_prev(prev.size()), n_curr(current.size());
		#pragma omp parallel for
		for (size_t p = 0; p < prev.size(); p++) {
			double _prev = 0, _curr = 0;
			for (size_t i = 0; i < current[p].size(); i++) {
				size_t commonDomains = mesh.getDOFsElement(p, i)->numberOfGlobalDomainsWithDOF(mesh.getDOFsElement(p, i)->DOFOffset(p, i));
				_curr += current[p][i] * current[p][i] / commonDomains;
				_prev += (current[p][i] - prev[p][i]) * (current[p][i] - prev[p][i]) / commonDomains;
			}
			n_prev[p] = _prev;
			n_curr[p] = _curr;
		}
		double localSum[2] = { std::accumulate(n_prev.begin(), n_prev.end(), 0), std::accumulate(n_curr.begin(), n_curr.end(), 0) };
		double globalSum[2];
		MPI_Allreduce(localSum, globalSum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return sqrt(globalSum[0]) / sqrt(globalSum[1]);
	};

	timeSolve.start();
	_linearSolver.Solve(_physics.f, solution);
	// reconstruct
	timeSolve.endWithBarrier();

	double error = 1;
	for (size_t i = 1; i < 100 && error > 0.01; i++) {
		#pragma omp parallel for
		for (size_t p = 0; p < _mesh.parts(); p++) {
			_prevSolution[p] = solution[p];
		}

		InequalityConstraints::removePositive(_constraints, solution, 1);
		// recompute GGt

		timeSolve.start();
		_linearSolver.Solve(_physics.f, solution);
		InequalityConstraints::reconstruct(_constraints);
		// reconstruct
		timeSolve.endWithBarrier();

		error = computeError(_mesh, _prevSolution, solution);
	}

	_timeStatistics.addEvent(timeSolve);

	if (_output.results) {
		_physics.saveMeshResults(_store, solution);
	}
}

template <class TPhysics, class TConfiguration>
void SemiSmoothNewtonMethod<TPhysics, TConfiguration>::finalize()
{
	_linearSolver.finilize();

	_timeStatistics.totalTime.endWithBarrier();
	_timeStatistics.printStatsMPI();
}

}
