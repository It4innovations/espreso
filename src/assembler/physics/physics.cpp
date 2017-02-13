
#include "physics.h"

#include "../../basis/logging/logging.h"
#include "../instance.h"

#include "../../mesh/elements/element.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"

#include "../constraints/equalityconstraints.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../basis/matrices/sparseVVPMatrix.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/matrices/sparseCSRMatrix.h"
#include "../../configuration/solver/espresooptions.h"


using namespace espreso;

Physics::Physics(Mesh *mesh, Instance *instance)
: _mesh(mesh), _instance(instance)
{

}


void Physics::assembleStiffnessMatrices(const Step &step)
{
	#pragma omp parallel for
	for  (size_t d = 0; d < _instance->domains; d++) {

		assembleStiffnessMatrix(step, d);

		switch (_instance->K[d].mtype) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			_instance->K[d].RemoveLower();
			break;
		case MatrixType::REAL_UNSYMMETRIC:
			break;
		}

		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

void Physics::assembleStiffnessMatrix(const Step &step, size_t domain)
{
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke, fe;
	std::vector<eslocal> DOFs;

	_K.resize(_instance->DOFs[domain], _instance->DOFs[domain]);
	_instance->f[domain].resize(_instance->DOFs[domain]);

	for (eslocal e = _mesh->getPartition()[domain]; e < _mesh->getPartition()[domain + 1]; e++) {
		processElement(step, _mesh->elements()[e], Ke, fe);
		fillDOFsIndices(_mesh->elements()[e], domain, DOFs);
		insertElementToDomain(_K, DOFs, Ke, fe, domain);
	}

	for (size_t i = 0; i < _mesh->faces().size(); i++) {
		if (_mesh->faces()[i]->inDomain(domain)) {
			processFace(step, _mesh->faces()[i], Ke, fe);
			fillDOFsIndices(_mesh->faces()[i], domain, DOFs);
			insertElementToDomain(_K, DOFs, Ke, fe, domain);
		}
	}

	for (size_t i = 0; i < _mesh->edges().size(); i++) {
		if (_mesh->edges()[i]->inDomain(domain)) {
			processEdge(step, _mesh->edges()[i], Ke, fe);
			fillDOFsIndices(_mesh->edges()[i], domain, DOFs);
			insertElementToDomain(_K, DOFs, Ke, fe, domain);
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	_instance->K[domain] = csrK;

	_instance->K[domain].mtype = getMatrixType(step, domain);
}

void Physics::subtractResidualForces(const Step &step)
{
	#pragma omp parallel for
	for  (size_t d = 0; d < _instance->domains; d++) {
		subtractResidualForces(step, d);
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

void Physics::subtractResidualForces(const Step &step, size_t domain)
{
	DenseMatrix Re;
	std::vector<eslocal> DOFs;

	for (eslocal e = _mesh->getPartition()[domain]; e < _mesh->getPartition()[domain + 1]; e++) {
		assembleResidualForces(step, _mesh->elements()[e], Re);
		fillDOFsIndices(_mesh->elements()[e], domain, DOFs);
		for (size_t i = 0; i < Re.rows(); i++) {
			_instance->f[domain][DOFs[i]] -= Re(i, 0);
		}
	}
}

void Physics::assembleStiffnessMatrix(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const
{
	size_t domain = e->domains().front();

	processElement(step, e, Ke, fe);

	for (size_t i = 0; i < e->filledFaces(); i++) {
		DenseMatrix Ki, fi;
		processFace(step, e->face(i), Ki, fi);
		Ke += Ki;
		fe += fi;
	}
	for (size_t i = 0; i < e->filledEdges(); i++) {
		DenseMatrix Ki, fi;
		processEdge(step, e->edge(i), Ki, fi);
		Ke += Ki;
		fe += fi;
	}
}

/**
 *
 * The method assumed that element matrix is composed in the following order:
 * x1, x2, x3, ..., y1, y2, y3, ..., z1, z2, z3,...
 *
 */
void Physics::fillDOFsIndices(const Element *e, eslocal domain, std::vector<eslocal> &DOFs) const
{
	DOFs.resize(e->nodes() * pointDOFs().size());
	for (size_t n = 0, i = 0; n < e->nodes(); n++) {
		for (size_t dof = 0; dof < pointDOFs().size(); dof++, i++) {
			DOFs[i] = _mesh->nodes()[e->node(n)]->DOFIndex(domain, dof);
		}
	}
}

void Physics::insertElementToDomain(SparseVVPMatrix<eslocal> &K, const std::vector<eslocal> &DOFs, const DenseMatrix &Ke, const DenseMatrix &fe, size_t domain)
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

void Physics::makeStiffnessMatricesRegular(REGULARIZATION regularization)
{
	#pragma omp parallel for
	for (size_t d = 0; d < _instance->domains; d++) {

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

double Physics::computeNormOfSolution() const
{
	double norm = 0, sum;
	double solution;
	for (size_t n = 0; n < _mesh->nodes().size(); n++) {
		for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
			size_t multiplicity = _mesh->nodes()[n]->numberOfGlobalDomainsWithDOF(dof);
			for (size_t d = 0; d < _mesh->nodes()[n]->domains().size(); d++) {
				size_t domain = _mesh->nodes()[n]->domains()[d];
				solution = _instance->primalSolution[domain][_mesh->nodes()[n]->DOFIndex(domain, dof)] / multiplicity;
				norm += solution * solution;
			}
		}
	}
	MPI_Allreduce(&norm, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(sum);
}

void Physics::assembleB1(const Step &step, bool withRedundantMultipliers, bool withScaling)
{
	EqualityConstraints::insertDirichletToB1(*_instance, _mesh->regions(), _mesh->nodes(), pointDOFs(), withRedundantMultipliers);
	EqualityConstraints::insertElementGluingToB1(*_instance, _mesh->neighbours(), _mesh->regions(), _mesh->nodes(), pointDOFs(), withRedundantMultipliers, withScaling);
}


