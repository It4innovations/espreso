
#include "linear.h"

namespace assembler {

template <class TInput>
void Linear<TInput>::init()
{
	this->_timeStatistics.totalTime.startWithBarrier();
	std::cout.precision(15);

	TimeEvent timeKasm("Create K and RHS");
	timeKasm.start();

	_K.resize(this->subdomains());
	_T.resize(this->subdomains());
	_M.resize(this->subdomains());
	_f.resize(this->subdomains());

	if (this->_verbose && esconfig::MPIrank == 0) {
		std::cout << "Assembling matrices : ";
	}
	cilk_for (size_t s = 0; s < this->subdomains(); s++) {
		//std::cout << s << " " ;
		// TODO: set dynamics
		KMf(s, false);
		T(s);


		if (this->_verbose && esconfig::MPIrank == 0) {
			std::cout << "." ;//<< s << " " ;
		}
	}
	if (this->_verbose && esconfig::MPIrank == 0) {
		std::cout << std::endl;
	}

	timeKasm.endWithBarrier();
	this->_timeStatistics.addEvent(timeKasm);

	TimeEvent timeLocalB("Create local B");
	timeLocalB.start();

	this->computeSubdomainGluing();

	timeLocalB.endWithBarrier();
	this->_timeStatistics.addEvent(timeLocalB);

	TimeEvent timeGlobalB("Create global B");
	timeGlobalB.start();

	std::vector<size_t> rows(this->subdomains());
	for (size_t s = 0; s < this->subdomains(); s++) {
		rows[s] = _K[s].rows;
	}

	this->computeClusterGluing(rows);

	timeGlobalB.endWithBarrier();
	this->_timeStatistics.addEvent(timeGlobalB);

	TimeEvent timeBforces("Fill right hand side");
	timeBforces.start();

	RHS();

	if (esconfig::info::printMatrices) {
		for (size_t s = 0; s < this->subdomains(); s++) {
			std::ofstream osK(eslog::Logging::prepareFile(s, "K").c_str());
			osK << _K[s];
			osK.close();
			std::ofstream osT(eslog::Logging::prepareFile(s, "T").c_str());
			osT << _T[s];
			osT.close();

			std::ofstream osF(eslog::Logging::prepareFile(s, "f").c_str());
			osF << _f[s];
			osF.close();

			std::ofstream osB0(eslog::Logging::prepareFile(s, "B0").c_str());
			osB0 << this->_B0[s];
			osB0.close();

			std::ofstream osB1(eslog::Logging::prepareFile(s, "B1").c_str());
			osB1 << this->_B1[s];
			osB1.close();
		}
	}


	timeBforces.endWithBarrier();
	this->_timeStatistics.addEvent(timeBforces);

	TimeEvent timeLSconv(string("Linear Solver - preprocessing"));
	timeLSconv.start();

	_lin_solver.DOFS_PER_NODE = this->DOFs();
	_lin_solver.setup(esconfig::MPIrank, esconfig::MPIsize, true);

	initSolver();

	timeLSconv.endWithBarrier();
	this->_timeStatistics.addEvent(timeLSconv);
}

template <class TInput>
void Linear<TInput>::pre_solve_update()
{

}

template <class TInput>
void Linear<TInput>::post_solve_update()
{
//	TimeEvent timeSaveVTK("Solver - Save VTK");
//	timeSaveVTK.start();
//
//	saveResult();
//
//	timeSaveVTK.endWithBarrier();
//	this->_timeStatistics.addEvent(timeSaveVTK);
}

template <class TInput>
void Linear<TInput>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeLSrun("Linear Solver - runtime");
	timeLSrun.start();

	_lin_solver.Solve(_f, solution);

	timeLSrun.endWithBarrier();
	this->_timeStatistics.addEvent(timeLSrun);
}

template <class TInput>
void Linear<TInput>::finalize()
{
	_lin_solver.finilize();

	this->_timeStatistics.totalTime.endWithBarrier();
	this->_timeStatistics.printStatsMPI();
}

template <class TInput>
void Linear<TInput>::T(size_t part)
{
	SparseIJVMatrix<eslocal> _T(1, 1);

	this->_T[part] = _T;
}


}
