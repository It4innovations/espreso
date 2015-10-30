
#include "linear.h"

namespace assembler {

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::init()
{
	this->_timeStatistics.SetName("Linear Elasticity Solver Overall Timing");
	this->_timeStatistics.totalTime.AddStartWithBarrier();
	std::cout.precision(15);

	TimeEvent timeKasm("Create K and RHS");
	timeKasm.AddStart();

	_K.resize(this->subdomains());
	_M.resize(this->subdomains());
	_f.resize(this->subdomains());
	for (size_t s = 0; s < this->subdomains(); s++) {
		std::cout << s << " " ;
		// TODO: set dynamics
		KMf(s, false);

		if (this->_verbose && this->_mesh.rank() == 0) {
			std::cout << s << " " ;
		}
	}
	if (this->_verbose && this->_mesh.rank() == 0) {
		std::cout << std::endl;
	}

	timeKasm.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeKasm);

	TimeEvent timeLocalB("Create local B");
	timeLocalB.AddStart();

	this->computeSubdomainGluing();

	timeLocalB.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeLocalB);

	TimeEvent timeGlobalB("Create global B");
	timeGlobalB.AddStart();

	std::vector<size_t> rows(this->subdomains());
	for (size_t s = 0; s < this->subdomains(); s++) {
		rows[s] = _K[s].rows;
	}
	this->computeClusterGluing(rows);

	timeGlobalB.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeGlobalB);

	TimeEvent timeBforces("Fill right hand side");
	timeBforces.AddStart();

	RHS();

	timeBforces.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeBforces);

	TimeEvent timeLSconv(string("Linear Solver - preprocessing"));
	timeLSconv.AddStart();

	_lin_solver.DOFS_PER_NODE = this->DOFs();
	_lin_solver.setup(this->_mesh.rank(), this->_mesh.size(), true);

	initSolver();

	timeLSconv.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeLSconv);
}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::pre_solve_update()
{

}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::post_solve_update()
{
	TimeEvent timeSaveVTK("Solver - Save VTK");
	timeSaveVTK.AddStart();

	saveResult();

	timeSaveVTK.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeSaveVTK);
}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::solve()
{
	TimeEvent timeLSrun("Linear Solver - runtime");
	timeLSrun.AddStart();

	_lin_solver.Solve(_f, _prim_solution);

	timeLSrun.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeLSrun);
}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::finalize()
{
	_lin_solver.finilize();

	this->_timeStatistics.totalTime.AddEndWithBarrier();
	this->_timeStatistics.PrintStatsMPI();
}


}
