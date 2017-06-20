
#include "../../../old_physics/precomputed/singular/assembler.h"

#include "../../../../basis/matrices/sparseVVPMatrix.h"
#include "../../../../basis/matrices/sparseCSRMatrix.h"
#include "../../../../solver/generic/SparseMatrix.h"

#include "../../../../mesh/structures/mesh.h"
#include "../../../../mesh/elements/element.h"
#include "../../../../mesh/settings/property.h"

#include "../../../constraints/equalityconstraints.h"

namespace espreso {

void SingularSystem::prepareMeshStructures()
{
	matrixSize = _apimesh.distributeDOFsToDomains(matrixSize);
	_apimesh.computeDOFsDOFsCounters();

	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		_apimesh.computeFacesSharedByDomains();
	}

	_constraints.initMatrices(matrixSize);
}

void SingularSystem::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	ESINFO(GLOBAL_ERROR) << "Implement assembleStiffnessMatrix";
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain, size_t scSize)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain, scSize);
}

void SingularSystem::makeStiffnessMatricesRegular()
{
	ESINFO(PROGRESS3) << "Make stiffness matrices regular.";
	#pragma omp parallel for
	for  (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		K[subdomain].RemoveLower();
		algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain, _solverConfiguration.SC_SIZE);
		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
}

void SingularSystem::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _apimesh.DOFs(), { Property::UNKNOWN });
	EqualityConstraints::insertElementGluingToB1(_constraints, _apimesh.DOFs(), { Property::UNKNOWN }, K);
}

void SingularSystem::assembleB0()
{
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		EqualityConstraints::insertKernelsToB0(_constraints, _apimesh.faces(), _apimesh.DOFs(), R1);
	}
}

void SingularSystem::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	eslocal nK = matrixSize[subdomain];
	_K.resize(nK, nK);

	const std::vector<eslocal> &partition = _apimesh.getPartition();
	const std::vector<Element*> &DOFs = _apimesh.DOFs();

	size_t dofIndex = 0;

	for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

		for (size_t dx = 0; dx < _apimesh.elements()[e]->DOFsIndices().size(); dx++) {
			size_t row = DOFs[_apimesh.elements()[e]->DOFsIndices()[dx]]->DOFIndex(subdomain, dofIndex);
			for (size_t dy = 0; dy < _apimesh.elements()[e]->DOFsIndices().size(); dy++) {
				size_t column = DOFs[_apimesh.elements()[e]->DOFsIndices()[dy]]->DOFIndex(subdomain, dofIndex);
				_K(row, column) = _apimesh.elements()[e]->stiffnessMatrix()[dx * _apimesh.elements()[e]->DOFsIndices().size() + dy];
			}
		}

	}

	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	f[subdomain].resize(matrixSize[subdomain], 0);
	for (size_t i = 0; i < DOFs.size(); i++) {
		if (DOFs[i]->inDomain(subdomain)) {
			f[subdomain][DOFs[i]->DOFIndex(subdomain, dofIndex)] = _rhs[i] / DOFs[i]->domains().size();
		}
	}
}
}
