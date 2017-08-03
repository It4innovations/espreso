
#include "precomputed.h"

#include "../step.h"
#include "../instance.h"
#include "../solution.h"

#include "../constraints/equalityconstraints.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/elementtypes.h"
#include "../../mesh/elements/element.h"

#include "../../basis/matrices/sparseVVPMatrix.h"
#include "../../basis/matrices/sparseCSRMatrix.h"
#include "../../solver/generic/SparseMatrix.h"

using namespace espreso;

Precomputed::Precomputed(Mesh *mesh, Instance *instance, MatrixType type, double *rhs, size_t rhsSize)
: Physics("API", mesh, instance), _mtype(type), _rhs(rhs), _rhsSize(rhsSize)
{

}

std::vector<size_t> Precomputed::solutionsIndicesToStore() const
{
	return { };
}

std::vector<std::pair<ElementType, Property> > Precomputed::propertiesToStore() const
{
	return { };
}

MatrixType Precomputed::getMatrixType(const Step &step, size_t domain) const
{
	return _mtype;
}

bool Precomputed::isMatrixTimeDependent(const Step &step) const
{
	return false;
}

bool Precomputed::isMatrixTemperatureDependent(const Step &step) const
{
	return false;
}

void Precomputed::prepareTotalFETI()
{
	_instance->domainDOFCount = dynamic_cast<APIMesh*>(_mesh)->distributeDOFsToDomains(_instance->domainDOFCount);
	dynamic_cast<APIMesh*>(_mesh)->computeDOFsDOFsCounters();
}

void Precomputed::prepareHybridTotalFETIWithKernels()
{
	prepareTotalFETI();
	dynamic_cast<APIMesh*>(_mesh)->computeFacesSharedByDomains();
}

void Precomputed::updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution)
{
	SparseVVPMatrix<eslocal> _K;

	if (matrices & Matrices::K) {
		_K.resize(_instance->domainDOFCount[domain], _instance->domainDOFCount[domain]);
	}
	if (matrices & Matrices::M) {
		ESINFO(GLOBAL_ERROR) << "PRECOMPUTED physics cannot assemble mass matrix M.";
	}
	if (matrices & Matrices::R) {
		ESINFO(GLOBAL_ERROR) << "PRECOMPUTED physics cannot assemble residual matrices.";
	}
	if (matrices & Matrices::f) {
		_instance->f[domain].clear();
		_instance->f[domain].resize(_instance->domainDOFCount[domain]);
	}

	const std::vector<eslocal> &partition = dynamic_cast<APIMesh*>(_mesh)->getPartition();
	const std::vector<Element*> &DOFs = dynamic_cast<APIMesh*>(_mesh)->DOFs();

	for (eslocal e = partition[domain]; e < partition[domain + 1]; e++) {

		for (size_t dx = 0; dx < dynamic_cast<APIMesh*>(_mesh)->elements()[e]->DOFsIndices().size(); dx++) {
			size_t row = DOFs[dynamic_cast<APIMesh*>(_mesh)->elements()[e]->DOFsIndices()[dx]]->DOFIndex(domain, 0);
			for (size_t dy = 0; dy < dynamic_cast<APIMesh*>(_mesh)->elements()[e]->DOFsIndices().size(); dy++) {
				size_t column = DOFs[dynamic_cast<APIMesh*>(_mesh)->elements()[e]->DOFsIndices()[dy]]->DOFIndex(domain, 0);
				_K(row, column) = dynamic_cast<APIMesh*>(_mesh)->elements()[e]->stiffnessMatrix()[dx * dynamic_cast<APIMesh*>(_mesh)->elements()[e]->DOFsIndices().size() + dy];
			}
		}
	}

	SparseCSRMatrix<eslocal> csrK = _K;
	_instance->K[domain] = csrK;
	_instance->K[domain].mtype = getMatrixType(step, domain);

	for (size_t i = 0; i < DOFs.size(); i++) {
		if (DOFs[i]->inDomain(domain)) {
			_instance->f[domain][DOFs[i]->DOFIndex(domain, 0)] = _rhs[i] / DOFs[i]->domains().size();
		}
	}
}

void Precomputed::assembleB1(const Step &step, bool withRedundantMultipliers, bool withScaling)
{
	EqualityConstraints::insertDirichletToB1(*_instance, step, _mesh->regions(), dynamic_cast<APIMesh*>(_mesh)->DOFs(), { Property::UNKNOWN }, { 0 }, withRedundantMultipliers);
	EqualityConstraints::insertElementGluingToB1(*_instance, step, _mesh->neighbours(), _mesh->regions(), dynamic_cast<APIMesh*>(_mesh)->DOFs(), { Property::UNKNOWN }, { 0 }, withRedundantMultipliers, withScaling);
}

void Precomputed::assembleB0FromKernels(const Step &step, const std::vector<SparseMatrix> &kernels)
{
	EqualityConstraints::insertKernelsGluingToB0(*_instance, _mesh->faces(), dynamic_cast<APIMesh*>(_mesh)->DOFs(), { 0 }, kernels, true);
}

void Precomputed::preprocessData(const Step &step)
{
	_instance->solutions.resize(SolutionIndex::SIZE, NULL);
	_instance->solutions[SolutionIndex::DECOMPOSED] = new Solution(*_mesh, "decomposed", ElementType::NODES, { Property::UNKNOWN }, _instance->primalSolution);
	_mergedSolution = std::vector<std::vector<double> > (1, std::vector<double>(dynamic_cast<APIMesh*>(_mesh)->DOFs().size(), 0));
	_instance->solutions[SolutionIndex::MERGED] = new Solution(*_mesh, "mergedSolution", ElementType::NODES, { Property::UNKNOWN }, _mergedSolution);
}

void Precomputed::processSolution(const Step &step)
{
	for (size_t i = 0; i < dynamic_cast<APIMesh*>(_mesh)->DOFs().size(); i++) {
		_mergedSolution[0][i] = 0;
		for (size_t d = 0; d < dynamic_cast<APIMesh*>(_mesh)->DOFs()[i]->domains().size(); d++) {
			eslocal domain = dynamic_cast<APIMesh*>(_mesh)->DOFs()[i]->domains()[d];
			eslocal dof = dynamic_cast<APIMesh*>(_mesh)->DOFs()[i]->DOFIndex(domain, 0);
			_mergedSolution[0][i] += _instance->solutions[SolutionIndex::DECOMPOSED]->data[domain][dof];
		}
		_mergedSolution[0][i] /= dynamic_cast<APIMesh*>(_mesh)->DOFs()[i]->domains().size();
	}
}


void Precomputed::prepareHybridTotalFETIWithCorners()
{
	ESINFO(ERROR) << "Cannot compute corners. Use HYBRID FETI with kernels.";
}

void Precomputed::assembleB0FromCorners(const Step &step)
{
	ESINFO(ERROR) << "Cannot compute corners. Use HYBRID FETI with kernels.";
}

void Precomputed::analyticRegularization(size_t domain)
{
	ESINFO(ERROR) << "Cannot compute analytic regularization of not PRECOMPUTED physics. Set REGULARIZATION = NULL_PIVOTS";
}

void Precomputed::processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	ESINFO(ERROR) << "ESPRESO internal error: cannot process element of precomputed physics";
}

void Precomputed::processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	ESINFO(ERROR) << "ESPRESO internal error: cannot process face of precomputed physics";
}

void Precomputed::processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	ESINFO(ERROR) << "ESPRESO internal error: cannot process edge of precomputed physics";
}

void Precomputed::processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	ESINFO(ERROR) << "ESPRESO internal error: cannot process node of precomputed physics";
}

