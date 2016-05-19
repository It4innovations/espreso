
#include "equalityconstraints.h"

namespace espreso {

template <class TInput>
EqualityConstraints<TInput>::EqualityConstraints(TInput &input): Assembler<TInput>(input)
{
	_B0.resize(this->subdomains());
	_B0subdomainsMap.resize(this->subdomains());

	_B1.resize(this->subdomains());
	_B1subdomainsMap.resize(this->subdomains());
	_B1duplicity.resize(this->subdomains());
	_B1c.resize(this->subdomains());

	for (size_t s = 0; s < this->subdomains(); s++) {
		_B0[s].rows = 0;
		_B0[s].nnz = 0;
		_B0[s].type = 'G';

		_B1[s].rows = 0;
		_B1[s].nnz = 0;
		_B1[s].type = 'G';
	}
}

static void initColumns(std::vector<SparseMatrix> &B0, std::vector<SparseMatrix> &B1, std::vector<size_t> &columns)
{
	for (size_t s = 0; s < columns.size(); s++) {
		B0[s].cols = columns[s];
		B1[s].cols = columns[s];
	}
}

static void postProcess(
		std::vector<SparseMatrix> &_B0,
		std::vector<SparseMatrix> &_B1,
		std::vector<std::vector<esglobal> > &_B0subdomainsMap,
		std::vector<std::vector<esglobal> > &_B1subdomainsMap,
		std::vector<std::vector<double> > &_B1duplicity,
		std::vector<std::vector<double> > &_B1c)
{
	cilk_for (size_t p = 0; p < _B1subdomainsMap.size(); p++) {
		// _B1subdomainsMap is the same as row indices
		_B1subdomainsMap[p] = _B1[p].I_row_indices;
		std::for_each(_B1subdomainsMap[p].begin(), _B1subdomainsMap[p].end(), [] (esglobal &v) { v -= IJVMatrixIndexing; });
		_B1c[p].resize(_B1[p].I_row_indices.size(), 0);
		std::sort(_B1subdomainsMap[p].begin(), _B1subdomainsMap[p].end());

		// because we sort subdomains map, we need to sort duplicity
		std::vector<double> origin(_B1duplicity[p].begin(), _B1duplicity[p].end());
		for (size_t i = 0; i < _B1[p].I_row_indices.size(); i++) {
			auto it = std::lower_bound(_B1subdomainsMap[p].begin(), _B1subdomainsMap[p].end(), _B1[p].I_row_indices[i] - 1);
			_B1duplicity[p][it - _B1subdomainsMap[p].begin()] = origin[i];
		}

		// _B0subdomainsMap is the same as row indices
		_B0subdomainsMap[p] = _B0[p].I_row_indices;
		std::for_each(_B0subdomainsMap[p].begin(), _B0subdomainsMap[p].end(), [] (esglobal &v) { v -= IJVMatrixIndexing; });
	}
}

template <>
void EqualityConstraints<API>::assembleConstraints(std::vector<size_t> columns)
{
	initColumns(_B0, _B1, columns);

	size_t lambdaCounter = 0;


	Dirichlet dir(
			*(this->_input.mesh), lambdaCounter,
			this->_input.dirichlet_size, this->_input.indexing,
			this->_input.dirichlet_indices, this->_input.dirichlet_values);

	lambdaCounter += dir.assemble(_B1, _B1clustersMap, _B1c);

	// TODO: duplicity is not important for dirichlet -> it is always 1
	cilk_for (size_t p = 0; p < this->_input.mesh->parts(); p++) {
		_B1duplicity[p].clear();
		_B1duplicity[p].resize(_B1[p].I_row_indices.size(), 1);
	}


	Gluing gluing(
			*(this->_input.mesh), lambdaCounter,
			this->_input.dirichlet_size, this->_input.indexing, this->_input.dirichlet_indices);

	std::vector<eslocal> corners;

	if (config::solver::FETI_METHOD == config::HYBRID_FETI && !config::solver::REDUNDANT_LAGRANGE) {
		// in this case B1 has to ignore corner nodes
		for (size_t i = 0; i < this->_input.mesh->subdomainBoundaries().size(); i++) {
			if (this->_input.mesh->subdomainBoundaries().isCorner(i)) {
				corners.push_back(i);
			}
		}
	}

	lambdaCounter += gluing.assembleB1(_B1, _B1clustersMap, _B1duplicity, corners);
	if (config::solver::FETI_METHOD == config::HYBRID_FETI) {
		gluing.assembleB0(_B0);
	}

	// TODO: get rid of the following
	postProcess(_B0, _B1, _B0subdomainsMap, _B1subdomainsMap, _B1duplicity, _B1c);
}

template <>
void EqualityConstraints<FEM>::assembleConstraints(std::vector<size_t> columns)
{
	initColumns(_B0, _B1, columns);

	// TODO: refactor DIRICHLET
	const std::map<eslocal, double> &dx = this->_input.mesh.coordinates().property(DIRICHLET_X).values();
	const std::map<eslocal, double> &dy = this->_input.mesh.coordinates().property(DIRICHLET_Y).values();
	const std::map<eslocal, double> &dz = this->_input.mesh.coordinates().property(DIRICHLET_Z).values();
	std::vector<eslocal> dirichlet;
	std::vector<double> dirichletValues;

	switch (this->_input.mesh.DOFs()) {
	case 1:
		dirichlet.reserve(dx.size());
		dirichletValues.reserve(dx.size());
		for (auto it = dx.begin(); it != dx.end(); ++it) {
			dirichlet.push_back(it->first);
			dirichletValues.push_back(it->second);
		}
		break;
	case 3:
		dirichlet.reserve(dx.size() + dy.size() + dz.size());
		dirichletValues.reserve(dx.size() + dy.size() + dz.size());
		for (auto it = dx.begin(); it != dx.end(); ++it) {
			dirichlet.push_back(3 * it->first);
			dirichletValues.push_back(0);
		}
		for (auto it = dy.begin(); it != dy.end(); ++it) {
			dirichlet.push_back(3 * it->first + 1);
			dirichletValues.push_back(0);
		}
		for (auto it = dz.begin(); it != dz.end(); ++it) {
			dirichlet.push_back(3 * it->first + 2);
			dirichletValues.push_back(0);
		}
		break;
	default:
		ESINFO(ERROR) << "Invalid number of DOFs";
	}

	const std::vector<BoundaryCondition*> &bc = this->_input.mesh.boundaryConditions();
	for (size_t i = 0; i < bc.size(); i++) {
		if (bc[i]->type() == ConditionType::DIRICHLET) {
			dirichlet.insert(dirichlet.end(), bc[i]->DOFs().begin(), bc[i]->DOFs().end());
			dirichletValues.insert(dirichletValues.end(), bc[i]->DOFs().size(), bc[i]->value());
		}
	}

	size_t lambdaCounter = 0;

	Dirichlet dir(this->_input.mesh, lambdaCounter, dirichlet, dirichletValues);
	lambdaCounter += dir.assemble(_B1, _B1clustersMap, _B1c);

	// TODO: duplicity is not important for dirichlet -> it is always 1
	cilk_for (size_t p = 0; p < this->_input.mesh.parts(); p++) {
		_B1duplicity[p].clear();
		_B1duplicity[p].resize(_B1[p].I_row_indices.size(), 1);
	}

	Gluing gluing(this->_input.mesh, lambdaCounter, dirichlet);

	std::vector<eslocal> corners;

	if (config::solver::FETI_METHOD == config::HYBRID_FETI && !config::solver::REDUNDANT_LAGRANGE) {
		// in this case B1 has to ignore corner nodes
		for (size_t i = 0; i < this->_input.mesh.subdomainBoundaries().size(); i++) {
			if (this->_input.mesh.subdomainBoundaries().isCorner(i)) {
				corners.push_back(i);
			}
		}
	}

	lambdaCounter += gluing.assembleB1(_B1, _B1clustersMap, _B1duplicity, corners);
	if (config::solver::FETI_METHOD == config::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::B0Type::CORNERS:
			gluing.assembleB0(_B0);
			break;
		case config::B0Type::KERNELS:
			gluing.assembleB0fromKernels(_B0);
			break;
		}
	}

	// TODO: get rid of the following
	postProcess(_B0, _B1, _B0subdomainsMap, _B1subdomainsMap, _B1duplicity, _B1c);
}

template <>
void EqualityConstraints<BEM>::assembleConstraints(std::vector<size_t> columns)
{
	initColumns(_B0, _B1, columns);

	// TODO: refactor DIRICHLET
	const std::map<eslocal, double> &dx = this->_input.mesh.coordinates().property(DIRICHLET_X).values();
	const std::map<eslocal, double> &dy = this->_input.mesh.coordinates().property(DIRICHLET_Y).values();
	const std::map<eslocal, double> &dz = this->_input.mesh.coordinates().property(DIRICHLET_Z).values();
	std::vector<eslocal> dirichlet;
	std::vector<double> dirichletValues;

	switch (this->_input.mesh.DOFs()) {
	case 1:
		dirichlet.reserve(dx.size());
		dirichletValues.reserve(dx.size());
		for (auto it = dx.begin(); it != dx.end(); ++it) {
			dirichlet.push_back(it->first);
			dirichletValues.push_back(it->second);
		}
		break;
	case 3:
		dirichlet.reserve(dx.size() + dy.size() + dz.size());
		dirichletValues.reserve(dx.size() + dy.size() + dz.size());
		for (auto it = dx.begin(); it != dx.end(); ++it) {
			dirichlet.push_back(3 * it->first);
			dirichletValues.push_back(0);
		}
		for (auto it = dy.begin(); it != dy.end(); ++it) {
			dirichlet.push_back(3 * it->first + 1);
			dirichletValues.push_back(0);
		}
		for (auto it = dz.begin(); it != dz.end(); ++it) {
			dirichlet.push_back(3 * it->first + 2);
			dirichletValues.push_back(0);
		}
		break;
	default:
		ESINFO(ERROR) << "Invalid number of DOFs";
	}

	size_t lambdaCounter = 0;

	Dirichlet dir(this->_input.surface, lambdaCounter, dirichlet, dirichletValues);
	lambdaCounter += dir.assemble(_B1, _B1clustersMap, _B1c);

	// TODO: duplicity is not important for dirichlet -> it is always 1
	cilk_for (size_t p = 0; p < this->_input.surface.parts(); p++) {
		_B1duplicity[p].clear();
		_B1duplicity[p].resize(_B1[p].I_row_indices.size(), 1);
	}

	Gluing gluing(this->_input.surface, lambdaCounter, dirichlet);

	std::vector<eslocal> corners;

	if (config::solver::FETI_METHOD == config::HYBRID_FETI && !config::solver::REDUNDANT_LAGRANGE) {
		// in this case B1 has to ignore corner nodes
		for (size_t i = 0; i < this->_input.mesh.subdomainBoundaries().size(); i++) {
			if (this->_input.mesh.subdomainBoundaries().isCorner(i)) {
				corners.push_back(i);
			}
		}
	}

	lambdaCounter += gluing.assembleB1(_B1, _B1clustersMap, _B1duplicity, corners);
	if (config::solver::FETI_METHOD == config::HYBRID_FETI) {
		gluing.assembleB0(_B0);
	}

	// TODO: get rid of the following
	postProcess(_B0, _B1, _B0subdomainsMap, _B1subdomainsMap, _B1duplicity, _B1c);
}


}
