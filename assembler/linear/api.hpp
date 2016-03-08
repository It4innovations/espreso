
#include "linear.h"

namespace assembler {

template <>
void Linear<API>::KMf(size_t part, bool dynamics)
{
	size_t DOFS = 3;
	SparseVVPMatrix<eslocal> _K;
	eslocal nK = this->_input.mesh->coordinates().localSize(part) * DOFS;
	_K.resize(nK, nK);

	const std::vector<eslocal> &parts = this->_input.mesh->getPartition();
	const std::vector<mesh::Element*> &elements = this->_input.mesh->getElements();
	const std::vector<std::vector<double> > &matrices = this->_input.mesh->getMatrices();

	size_t row, column;

	for (size_t e = parts[part]; e < parts[part + 1]; e++) {
		for (size_t i = 0; i < elements[e]->size() * DOFS; i++) {
			row = DOFS * elements[e]->node(i / DOFS) + (i % DOFS);
			for (size_t j = 0; j < elements[e]->size() * DOFS; j++) {
				column = DOFS * elements[e]->node(j / DOFS) + (j % DOFS);
				_K(row, column) = matrices[e][i * elements[e]->size() * DOFS + j];
			}
		}
	}

	SparseCSRMatrix<eslocal> csrK = _K;
	this->_K[part] = csrK;
}

template <>
void Linear<API>::RHS()
{
	size_t DOFS = 3;

	for (size_t p = 0; p < this->_input.mesh->parts(); p++) {
		const std::vector<eslocal> &l2c = this->_input.mesh->coordinates().localToCluster(p);
		const mesh::Boundaries &boundaries = this->_input.mesh->subdomainBoundaries();
		_f[p].resize(this->_input.mesh->coordinates().localSize(p) * DOFS, 0);
		for (size_t i = 0; i < l2c.size() * DOFS; i++) {
			_f[p][i] = this->_input.rhs[DOFS * l2c[i / DOFS] + i % DOFS] / boundaries[l2c[i / DOFS]].size();
		}
	}
}

template <>
void Linear<API>::initSolver()
{
	// TODO: solve it better

	std::vector<std::vector<eslocal> > fixPoints;

	_lin_solver.init(
		*_input.mesh,
		_K,
		_T,
		_B1,
		_B0,
		_B1subdomainsMap,
		_B0subdomainsMap,
		_B1clustersMap,
		_B1duplicity,
		_f,
		_B1c,
		fixPoints,
		_input.mesh->neighbours()
	);
}

template <>
void Linear<API>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeLSrun("Linear Solver - runtime");
	timeLSrun.start();

	std::vector<std::vector<double> > tmpSolution(_f.size());
	for (size_t i = 0; i <_f.size(); i++) {
		tmpSolution[i].resize(_f[i].size());
	}

	_lin_solver.Solve(_f, tmpSolution);

	size_t DOFs = 3;

	std::for_each(solution[0].begin(), solution[0].end(), [] (double &v) { v = 0; });

	const mesh::Boundaries &boundaries = this->_input.mesh->subdomainBoundaries();
	for (size_t p = 0; p < this->_input.mesh->parts(); p++) {
		const std::vector<eslocal> &l2c = this->_input.mesh->coordinates().localToCluster(p);
		for (size_t i = 0; i < l2c.size(); i++) {
			for (size_t d = 0; d < DOFs; d++) {
				solution[0][DOFs * l2c[i] + d] += tmpSolution[p][DOFs * i + d] / boundaries[l2c[i]].size();
			}
		}
	}


	timeLSrun.endWithBarrier();
	this->_timeStatistics.addEvent(timeLSrun);
}

}
