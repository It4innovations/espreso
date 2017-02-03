
#include "physics.h"

#include "../../basis/logging/logging.h"
#include "../../config/solverespresooptions.h"

#include "../instance.h"

#include "../../mesh/elements/element.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"

#include "../constraints/equalityconstraints.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../basis/matrices/sparseVVPMatrix.h"
#include "../../basis/matrices/denseMatrix.h"


using namespace espreso;

NewPhysics::NewPhysics(Mesh *mesh, NewInstance *instance)
: _mesh(mesh), _instance(instance)
{

}


void NewPhysics::assembleStiffnessMatrices()
{
	#pragma omp parallel for
	for  (size_t p = 0; p < _instance->domains; p++) {
		assembleStiffnessMatrix(p);
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

/**
 *
 * The method assumed that element matrix is composed in the following order:
 * x1, x2, x3, ..., y1, y2, y3, ..., z1, z2, z3,...
 *
 */
void NewPhysics::fillDOFsIndices(const Element *e, eslocal domain, std::vector<eslocal> &DOFs)
{
	DOFs.resize(e->nodes() * pointDOFs().size());
	for (size_t n = 0, i = 0; n < e->nodes(); n++) {
		for (size_t dof = 0; dof < pointDOFs().size(); dof++, i++) {
			DOFs[i] = _mesh->nodes()[e->node(n)]->DOFIndex(domain, dof);
		}
	}
}

void NewPhysics::insertElementToDomain(SparseVVPMatrix<eslocal> &K, const std::vector<eslocal> &DOFs, const DenseMatrix &Ke, const DenseMatrix &fe, size_t domain)
{
	if (Ke.rows() == DOFs.size() && Ke.columns() == DOFs.size()) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			for (size_t c = 0; c < DOFs.size(); c++) {
				K(DOFs[r], DOFs[c]) = Ke(r, c);
			}
		}
	} else {
		if (Ke.rows() != 0 || Ke.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while inserting element matrix to domain matrix";
		}
	}

	if (fe.rows() == DOFs.size() && fe.columns() == 1) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			_instance->f[domain][DOFs[r]] += fe(r, 0);
		}
	} else {
		if (fe.rows() != 0 || fe.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while inserting element RHS to domain RHS";
		}
	}
}

void NewPhysics::makeStiffnessMatricesRegular(REGULARIZATION regularization)
{
	#pragma omp parallel for
	for (size_t d = 0; d < _instance->domains; d++) {
		switch (_instance->K[d].mtype) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			_instance->K[d].RemoveLower();
			break;
		case MatrixType::REAL_UNSYMMETRIC:
			break;
		}

		switch (regularization) {

		case REGULARIZATION::FIX_POINTS:
			analyticRegularization(d);
			_instance->RegMat[d].RemoveLower();
			_instance->K[d].MatAddInPlace(_instance->RegMat[d], 'N', 1);
			_instance->RegMat[d].ConvertToCOO(1);
			break;

		case REGULARIZATION::NULL_PIVOTS:
			switch (_instance->K[d].mtype) {
				double norm;
				eslocal defect;

			case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
				_instance->K[d].get_kernel_from_K(_instance->K[d], _instance->RegMat[d], _instance->R1[d], norm, defect, d);
				break;

			case MatrixType::REAL_UNSYMMETRIC:
				_instance->K[d].get_kernels_from_nonsym_K(_instance->K[d], _instance->RegMat[d], _instance->R1[d], _instance->R2[d], norm, defect, d);
				break;

			default:
				ESINFO(ERROR) << "Unknown matrix type for regularization.";
			}
			break;
		}
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

void NewPhysics::assembleB1(bool withRedundantMultipliers, bool withScaling)
{
	EqualityConstraints::insertDirichletToB1(*_instance, _mesh->regions(), _mesh->nodes(), pointDOFs(), withRedundantMultipliers);
	EqualityConstraints::insertElementGluingToB1(*_instance, _mesh->neighbours(), _mesh->regions(), _mesh->nodes(), pointDOFs(), withRedundantMultipliers, withScaling);
}


