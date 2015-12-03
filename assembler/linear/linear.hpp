
#include "linear.h"

namespace assembler {

template <class TInput>
void Linear<TInput>::init()
{
	this->_timeStatistics.SetName("Linear Elasticity Solver Overall Timing");
	this->_timeStatistics.totalTime.AddStartWithBarrier();
	std::cout.precision(15);

	TimeEvent timeKasm("Create K and RHS");
	timeKasm.AddStart();

	_K.resize(this->subdomains());
	_M.resize(this->subdomains());
	_f.resize(this->subdomains());
	cilk_for (size_t s = 0; s < this->subdomains(); s++) {
		std::cout << s << " " ;
		// TODO: set dynamics
		KMf(s, false);

		if (this->_verbose && this->rank() == 0) {
			std::cout << s << " " ;
		}
	}
	if (this->_verbose && this->rank() == 0) {
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
	_lin_solver.setup(this->rank(), this->size(), true);

	initSolver();

	timeLSconv.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeLSconv);
}

template <class TInput>
void Linear<TInput>::pre_solve_update()
{

}

template <class TInput>
void Linear<TInput>::post_solve_update()
{
//	TimeEvent timeSaveVTK("Solver - Save VTK");
//	timeSaveVTK.AddStart();
//
//	saveResult();
//
//	timeSaveVTK.AddEndWithBarrier();
//	this->_timeStatistics.AddEvent(timeSaveVTK);
}

template <class TInput>
void Linear<TInput>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeLSrun("Linear Solver - runtime");
	timeLSrun.AddStart();

	_lin_solver.Solve(_f, solution);

	timeLSrun.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeLSrun);
}

template <class TInput>
void Linear<TInput>::finalize()
{
	_lin_solver.finilize();

	this->_timeStatistics.totalTime.AddEndWithBarrier();
	this->_timeStatistics.PrintStatsMPI();
}

template <>
void Linear<API>::fillAPIHolder(APIHolder *holder)
{

}

template <class TInput>
void Linear<TInput>::fillAPIHolder(APIHolder *holder)
{
	eslocal indexing = 1;

	init();

	SparseVVPMatrix<eslocal> vvp(_K[0].rows, _K[0].cols);
	for (size_t i = 0; i < _K[0].rows; i++) {
		for (size_t j = _K[0].CSR_I_row_indices[i]; j < _K[0].CSR_I_row_indices[i + 1]; j++) {
			vvp(i, _K[0].CSR_J_col_indices[j - 1] - 1) = _K[0].CSR_V_values[j - 1];
		}
	}
	holder->K = new SparseCSRMatrix<eslocal>(vvp);

	holder->rhs = new ESPRESOStructDoubleVector();
	holder->rhs->size = _f[0].size();
	holder->rhs->values = new double[_f[0].size()];
	for (size_t i = 0; i < _f[0].size(); i++) {
		holder->rhs->values[i] = _f[0][i];
	}

	std::vector<std::pair<esglobal, double> > dir;
	const mesh::Coordinates &coo = this->_input.mesh.coordinates();
	const std::map<eslocal, double> &dx = coo.property(mesh::DIRICHLET_X).values();
	const std::map<eslocal, double> &dy = coo.property(mesh::DIRICHLET_X).values();
	const std::map<eslocal, double> &dz = coo.property(mesh::DIRICHLET_X).values();
	for (size_t i = 0; i < coo.size(); i++) {
		if (dx.find(i) != dx.end()) {
			dir.push_back(std::pair<esglobal, double>(3 * coo.globalIndex(i) + indexing, dx.find(i)->second));
		}
		if (dy.find(i) != dy.end()) {
			dir.push_back(std::pair<esglobal, double>(3 * coo.globalIndex(i) + 1 + indexing, dy.find(i)->second));
		}
		if (dz.find(i) != dz.end()) {
			dir.push_back(std::pair<esglobal, double>(3 * coo.globalIndex(i) + 2 + indexing, dz.find(i)->second));
		}
	}
	holder->dirichlet = new ESPRESOStructMap();
	holder->dirichlet->size = dir.size();
	holder->dirichlet->indices = new eslocal[dir.size()];
	holder->dirichlet->values = new double[dir.size()];
	for (size_t i = 0; i < dir.size(); i++) {
		holder->dirichlet->indices[i] = dir[i].first;
		holder->dirichlet->values[i] = dir[i].second;
	}

	holder->l2g = new ESPRESOStructIntVector();
	holder->l2g->size = coo.size() * 3;
	holder->l2g->values = new eslocal[coo.size() * 3];
	for (size_t i = 0; i < coo.size(); i++) {
		holder->l2g->values[3 * i] = 3 * coo.globalIndex(i) + indexing;
		holder->l2g->values[3 * i + 1] = 3 * coo.globalIndex(i) + 1 + indexing;
		holder->l2g->values[3 * i + 2] = 3 * coo.globalIndex(i) + 2 + indexing;
	}

	holder->neighbourRanks = new ESPRESOStructIntVector();
	holder->neighbourRanks->size = this->_neighClusters.size();
	holder->neighbourRanks->values = new eslocal[this->_neighClusters.size()];
	for (size_t i = 0; i < this->_neighClusters.size(); i++) {
		holder->neighbourRanks->values[i] = this->_neighClusters[i] + indexing;
	}

	holder->indexing = indexing;
}


}
