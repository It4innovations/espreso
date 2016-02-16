
#include "linear.h"

namespace assembler {

template <>
void Linear<API>::KMf(size_t part, bool dynamics)
{
	if (esconfig::mesh::subdomains == 1) {
		_K[part] = *(this->_input.K);
		return;
	}
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
	if (esconfig::mesh::subdomains == 1) {
		_f[0] = std::vector<double>(this->_input.rhs, this->_input.rhs + this->_input.size);
		return;
	}
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
	_lin_solver.init(
		_K,
		_globalB,
		_localB,
		_lambda_map_sub_B1,
		_lambda_map_sub_B0,
		_lambda_map_sub_clst,
		_B1_duplicity,
		_f,
		_vec_c,
		_neighClusters
	);
}

}
