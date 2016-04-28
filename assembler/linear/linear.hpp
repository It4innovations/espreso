
#include "linear.h"

namespace espreso {

template <class TInput>
void Linear<TInput>::init()
{
	this->_timeStatistics.totalTime.startWithBarrier();

	TimeEvent timeKasm("Create K and RHS");
	timeKasm.start();

	_K.resize(this->subdomains());
	_T.resize(this->subdomains());
	_M.resize(this->subdomains());
	_f.resize(this->subdomains());

	ESINFO(PROGRESS2) << "Assemble matrices K, M, T and right hand side";
	cilk_for (size_t s = 0; s < this->subdomains(); s++) {
		KMf(s, config::assembler::timeSteps > 1);
		if (config::assembler::timeSteps > 1) {
			_K[s].MatAddInPlace(_M[s], 'N', timeConstant());
		}
		T(s);
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);

	timeKasm.endWithBarrier();
	this->_timeStatistics.addEvent(timeKasm);

	std::vector<size_t> rows(this->subdomains());
	for (size_t s = 0; s < this->subdomains(); s++) {
		rows[s] = _K[s].rows;
	}

	TimeEvent timeBforces("Fill right hand side");
	timeBforces.start();

	RHS();

	TimeEvent timeParallelG("Gluing");
	timeParallelG.startWithBarrier();

	ESINFO(PROGRESS2) << "Assemble equality constraints";
	this->assembleConstraints(rows);

	timeParallelG.end();
	this->_timeStatistics.addEvent(timeParallelG);

	if (config::info::printMatrices) {
		for (size_t s = 0; s < this->subdomains(); s++) {
			std::ofstream osK(Logging::prepareFile(s, "K").c_str());
			osK << _K[s];
			osK.close();
			std::ofstream osT(Logging::prepareFile(s, "T").c_str());
			osT << _T[s];
			osT.close();

			std::ofstream osF(Logging::prepareFile(s, "f").c_str());
			osF << _f[s];
			osF.close();

			std::ofstream osB0(Logging::prepareFile(s, "B0").c_str());
			osB0 << this->_B0[s];
			osB0.close();

			std::ofstream osB1(Logging::prepareFile(s, "B1").c_str());
			osB1 << this->_B1[s];
			osB1.close();
		}
	}

	timeBforces.endWithBarrier();
	this->_timeStatistics.addEvent(timeBforces);

	TimeEvent timeLSconv(string("Linear Solver - preprocessing"));
	timeLSconv.start();

	_lin_solver.DOFS_PER_NODE = this->DOFs();
	_lin_solver.setup(config::env::MPIrank, config::env::MPIsize, config::assembler::timeSteps == 1);

	initSolver();

	timeLSconv.endWithBarrier();
	this->_timeStatistics.addEvent(timeLSconv);
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
