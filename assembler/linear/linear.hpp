
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

template <>
void Linear<API>::fillAPIHolder(APIHolder *holder)
{

}

template <class TInput>
void Linear<TInput>::fillAPIHolder(APIHolder *holder)
{
	//TODO: dissabled because of 64bit int
//	eslocal indexing = 0;
//
//	init();
//
//	SparseVVPMatrix<eslocal> vvp(_K[0].rows, _K[0].cols);
//	for (size_t i = 0; i < _K[0].rows; i++) {
//		for (size_t j = _K[0].CSR_I_row_indices[i]; j < _K[0].CSR_I_row_indices[i + 1]; j++) {
//			vvp(i, _K[0].CSR_J_col_indices[j - 1] - 1) = _K[0].CSR_V_values[j - 1];
//		}
//	}
//	holder->K = new SparseCSRMatrix<eslocal>(vvp);
//
//	holder->rhs = new std::vector<double>(_f[0]);
//
//	std::map<eslocal, double> dirichlet;
//	const mesh::Coordinates &coo = this->_input.mesh.coordinates();
//	const std::map<eslocal, double> &dx = coo.property(mesh::DIRICHLET_X).values();
//	const std::map<eslocal, double> &dy = coo.property(mesh::DIRICHLET_X).values();
//	const std::map<eslocal, double> &dz = coo.property(mesh::DIRICHLET_X).values();
//	std::map<eslocal, double>::const_iterator it;
//	for (it = dx.begin(); it != dx.end(); ++it) {
//		dirichlet[3 * coo.globalIndex(it->first) + indexing] = it->second;
//	}
//	for (it = dy.begin(); it != dy.end(); ++it) {
//		dirichlet[3 * coo.globalIndex(it->first) + 1 + indexing] = it->second;
//	}
//	for (it = dz.begin(); it != dz.end(); ++it) {
//		dirichlet[3 * coo.globalIndex(it->first) + 2 + indexing] = it->second;
//	}
//
//	holder->dirichlet_indices = new std::vector<eslocal>(dirichlet.size());
//	holder->dirichlet_values = new std::vector<double>(dirichlet.size());
//	eslocal index = 0;
//	for (it = dirichlet.begin(); it != dirichlet.end(); ++it, index++) {
//		(*holder->dirichlet_indices)[index] = it->first;
//		(*holder->dirichlet_values)[index] = it->second;
//	}
//
//	holder->l2g = new std::vector<eslocal>(coo.clusterSize() * 3);
//	for (size_t i = 0; i < coo.clusterSize(); i++) {
//		(*holder->l2g)[3 * i] = 3 * coo.globalIndex(i) + indexing;
//		(*holder->l2g)[3 * i + 1] = 3 * coo.globalIndex(i) + 1 + indexing;
//		(*holder->l2g)[3 * i + 2] = 3 * coo.globalIndex(i) + 2 + indexing;
//	}
//
//	holder->neighbourRanks = new std::vector<eslocal>(this->_neighClusters);
//
//	holder->indexing = indexing;
}


}
