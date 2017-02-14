
#include "physics.h"

#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"
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

void Physics::assembleResidualForces(const Step &step)
{
	#pragma omp parallel for
	for  (size_t d = 0; d < _instance->domains; d++) {
		assembleResidualForces(step, d);
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

void Physics::assembleResidualForces(const Step &step, size_t domain)
{
	DenseMatrix Re;
	std::vector<eslocal> DOFs;

	_instance->R[domain].resize(_instance->DOFs[domain]);

	for (eslocal e = _mesh->getPartition()[domain]; e < _mesh->getPartition()[domain + 1]; e++) {
		assembleResidualForces(step, _mesh->elements()[e], Re);
		fillDOFsIndices(_mesh->elements()[e], domain, DOFs);
		for (size_t i = 0; i < Re.rows(); i++) {
			_instance->R[domain][DOFs[i]] += Re(i, 0);
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
				_instance->K[d].get_kernel_from_K(_instance->K[d], _instance->RegMat[d], _instance->N1[d], norm, defect, d);
				break;

			case MatrixType::REAL_UNSYMMETRIC:
				_instance->K[d].get_kernels_from_nonsym_K(_instance->K[d], _instance->RegMat[d], _instance->N1[d], _instance->N2[d], norm, defect, d);
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

double Physics::sumSquares(const std::vector<std::vector<double> > &data, SumOperation operation, SumRestriction restriction, size_t loadStep) const
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh->neighbours().begin(), _mesh->neighbours().end(), neighbour) - _mesh->neighbours().begin();
	};

	double csum = 0, gsum;

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _mesh->nodes().size());

	std::vector<std::vector<std::vector<double> > > sBuffer(threads, std::vector<std::vector<double> >(_mesh->neighbours().size()));
	std::vector<std::vector<double> > rBuffer(_mesh->neighbours().size());
	std::vector<std::vector<size_t> > rBufferSize(threads, std::vector<size_t>(_mesh->neighbours().size()));
	std::vector<std::vector<eslocal> > incomplete(threads);
	std::vector<std::vector<double> > incompleteData(threads);

	std::vector<std::vector<Region*> > allowed(pointDOFs().size());
	if (restriction == SumRestriction::DIRICHLET) {
		allowed = _mesh->getRegionsWithProperties(loadStep, pointDOFs());
	}


	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		double tSum = 0;
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {

			if (_mesh->nodes()[n]->clusters().size() > 1) {
				if (_mesh->nodes()[n]->clusters().front() == environment->MPIrank) {
					for (auto c = _mesh->nodes()[n]->clusters().begin() + 1; c != _mesh->nodes()[n]->clusters().end(); ++c) {
						rBufferSize[t][n2i(*c)] += pointDOFs().size();
					}
					incomplete[t].push_back(n);
					incompleteData[t].resize(incompleteData.size() + pointDOFs().size());
					for (auto d = _mesh->nodes()[n]->domains().begin(); d != _mesh->nodes()[n]->domains().end(); ++d) {
						for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
							if (restriction == SumRestriction::NONE || _mesh->commonRegion(allowed[dof], _mesh->nodes()[n]->regions())) {
								eslocal index = _mesh->nodes()[n]->DOFIndex(*d, dof);
								if (index >= 0) {
									incompleteData[t][incompleteData.size() - pointDOFs().size() + dof] += data[*d][index];
								}
							}
						}
					}
				} else {
					eslocal cluster = _mesh->nodes()[n]->clusters().front();
					sBuffer[t][n2i(cluster)].resize(sBuffer[t][n2i(cluster)].size() + pointDOFs().size());
					for (auto d = _mesh->nodes()[n]->domains().begin(); d != _mesh->nodes()[n]->domains().end(); ++d) {
						for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
							if (restriction == SumRestriction::NONE || _mesh->commonRegion(allowed[dof], _mesh->nodes()[n]->regions())) {
								eslocal index = _mesh->nodes()[n]->DOFIndex(*d, dof);
								if (index >= 0) {
									sBuffer[t][n2i(cluster)][sBuffer[t][n2i(cluster)].size() - pointDOFs().size() + dof] += data[*d][index];
								}
							}
						}
					}
				}
			} else {
				for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
					double sum = 0;
					for (auto d = _mesh->nodes()[n]->domains().begin(); d != _mesh->nodes()[n]->domains().end(); ++d) {
						if (restriction == SumRestriction::NONE || _mesh->commonRegion(allowed[dof], _mesh->nodes()[n]->regions())) {
							eslocal index = _mesh->nodes()[n]->DOFIndex(*d, dof);
							if (index >= 0) {
								sum += data[*d][index];
							}
						}
					}
					switch (operation) {
					case SumOperation::AVERAGE:
						sum /= _mesh->nodes()[n]->numberOfGlobalDomainsWithDOF(dof);
					case SumOperation::SUM:
						tSum += sum * sum;
						break;
					default:
						ESINFO(GLOBAL_ERROR) << "Implement new SumOperation.";
					}
				}
			}

		}
		csum += tSum;
	}

	for (size_t n = 0; n < _mesh->neighbours().size(); n++) {
		rBuffer[n].resize(rBuffer[n].size() + rBufferSize[0][n]);
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
			incomplete[0].insert(incomplete[0].end(), incomplete[t].begin(), incomplete[t].end());
			incompleteData[0].insert(incompleteData[0].end(), incompleteData[t].begin(), incompleteData[t].end());
			rBuffer[n].resize(rBuffer[n].size() + rBufferSize[t][n]);
		}
	}

	if (!Communication::receiveUpperKnownSize(sBuffer[0], rBuffer, _mesh->neighbours())) {
		ESINFO(ERROR) << "problem while exchange sum of squares.";
	}

	distribution = Esutils::getDistribution(threads, incomplete[0].size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		double tSum = 0;
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			size_t n = incomplete[0][i];

			for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
				double sum = incompleteData[0][n * pointDOFs().size() + dof];
				for (auto c = _mesh->nodes()[n]->clusters().begin() + 1; c != _mesh->nodes()[n]->clusters().end(); ++c) {
					sum += rBuffer[n2i(*c)][_mesh->nodes()[n]->clusterOffset(*c) * pointDOFs().size() + dof];
				}
				switch (operation) {
				case SumOperation::AVERAGE:
					sum /= _mesh->nodes()[n]->numberOfGlobalDomainsWithDOF(dof);
				case SumOperation::SUM:
					tSum += sum * sum;
					break;
				default:
					ESINFO(GLOBAL_ERROR) << "Implement new SumOperation.";
				}
			}

		}
		csum += tSum;
	}

	MPI_Allreduce(&csum, &gsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return gsum;
}

void Physics::assembleB1(const Step &step, bool withRedundantMultipliers, bool withScaling)
{
	EqualityConstraints::insertDirichletToB1(*_instance, _mesh->regions(), _mesh->nodes(), pointDOFs(), withRedundantMultipliers);
	EqualityConstraints::insertElementGluingToB1(*_instance, _mesh->neighbours(), _mesh->regions(), _mesh->nodes(), pointDOFs(), withRedundantMultipliers, withScaling);
}


