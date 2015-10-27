
#include "linear.h"

namespace physics {

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::init()
{
	std::cout << "0\n";
	this->_timeStatistics.SetName("Linear Elasticity Solver Overall Timing");
	this->_timeStatistics.totalTime.AddStartWithBarrier();
	std::cout.precision(15);

	TimeEvent timeKasm("Create K and RHS");
	timeKasm.AddStart();

	_K.resize(subdomains());
	_M.resize(subdomains());
	_f.resize(subdomains());
	for (size_t s = 0; s < subdomains(); s++) {
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
	std::cout << "1\n";

	TimeEvent timeLocalB("Create local B");
	timeLocalB.AddStart();

	_localB.resize(subdomains());
	_globalB.resize(subdomains());
	_B0.resize(subdomains());
	_B1.resize(subdomains());

	_lambda_map_sub_B1.resize(subdomains());
	_lambda_map_sub_B0.resize(subdomains());
	_B1_duplicity.resize(subdomains());
	_vec_c.resize(subdomains());
	localB();

	timeLocalB.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeLocalB);
	std::cout << "2\n";

	TimeEvent timeGlobalB("Create global B");
	timeGlobalB.AddStart();


	globalB();
	for (int i = 0; i < _vec_c.size(); i++) {
		 _vec_c[i].resize(_B1_duplicity[i].size(), 0.0);
	}

	for (size_t p = 0; p < this->_mesh.parts(); p++) {
		_localB[p] = _B0[p];
		_globalB[p] = _B1[p];
	}

	timeGlobalB.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeGlobalB);
	std::cout << "3\n";


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
	std::cout << "4\n";
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
