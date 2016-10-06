
#include "assembler.h"

using namespace espreso;

void UniformSymmetric3DOFs::prepareMeshStructures()
{
	matrixSize = _apimesh.distributeDOFsToDomains(matrixSize);
	_apimesh.computeDOFsDOFsCounters();

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		ESINFO(GLOBAL_ERROR) << "Implement Hybrid FETI for API";
	}
}

void UniformSymmetric3DOFs::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	ESINFO(GLOBAL_ERROR) << "Implement assembleStiffnessMatrix";
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void UniformSymmetric3DOFs::makeStiffnessMatricesRegular()
{
	cilk_for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		K[subdomain].RemoveLower();
		algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
	}
}

void UniformSymmetric3DOFs::assembleGluingMatrices()
{
	_constraints.initMatrices(matrixSize);

	_constraints.insertDirichletToB1(_apimesh.DOFs(), { Property::UNKNOWN });
	_constraints.insertElementGluingToB1(_apimesh.DOFs(), { Property::UNKNOWN });
}

void UniformSymmetric3DOFs::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	eslocal nK = matrixSize[subdomain];
	_K.resize(nK, nK);

	const std::vector<eslocal> &partition = _apimesh.getPartition();
	const std::vector<Element*> &DOFs = _apimesh.DOFs();

	size_t dofIndex = 0;

	for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

		for (size_t dx = 0; dx < _apimesh.eDOFs(e).size(); dx++) {
			size_t row = DOFs[_apimesh.eDOFs(e)[dx]]->DOFIndex(subdomain, dofIndex);
			for (size_t dy = 0; dy < _apimesh.eDOFs(e).size(); dy++) {
				size_t column = DOFs[_apimesh.eDOFs(e)[dy]]->DOFIndex(subdomain, dofIndex);
				_K(row, column) = _apimesh.eMatrix(e)[dx * _apimesh.eDOFs(e).size() + dy];
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
