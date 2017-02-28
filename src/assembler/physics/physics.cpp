
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

Physics::Physics(const std::string &name, Mesh *mesh, Instance *instance)
: _name(name), _mesh(mesh), _instance(instance)
{

}

void Physics::assembleMatrix(const Step &step, Matrices matrix)
{
	updateMatrix(step, matrix, {});
}

void Physics::updateMatrix(const Step &step, Matrices matrix, const std::vector<Solution*> &solution)
{
	#pragma omp parallel for
	for  (size_t d = 0; d < _instance->domains; d++) {

		updateMatrix(step, matrix, d, solution);

		switch (_instance->K[d].mtype) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			_instance->K[d].RemoveLower();
			_instance->M[d].RemoveLower();
			break;
		case MatrixType::REAL_UNSYMMETRIC:
			break;
		}

		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
}

void Physics::assembleMatrix(const Step &step, Matrices matrices, size_t domain)
{
	updateMatrix(step, matrices, domain, {});
}

void Physics::updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution)
{
	SparseVVPMatrix<eslocal> _K, _M;
	DenseMatrix Ke, Me, Re, fe;
	std::vector<eslocal> DOFs;

	if (matrices & Matrices::K) {
		_K.resize(_instance->DOFs[domain], _instance->DOFs[domain]);
	}
	if (matrices & Matrices::M) {
		_M.resize(_instance->DOFs[domain], _instance->DOFs[domain]);
	}
	if (matrices & Matrices::R) {
		_instance->R[domain].clear();
		_instance->R[domain].resize(_instance->DOFs[domain]);
	}
	if (matrices & Matrices::f) {
		_instance->f[domain].clear();
		_instance->f[domain].resize(_instance->DOFs[domain]);
	}

	for (eslocal e = _mesh->getPartition()[domain]; e < _mesh->getPartition()[domain + 1]; e++) {
		processElement(step, matrices, _mesh->elements()[e], Ke, Me, Re, fe, solution);
		fillDOFsIndices(_mesh->elements()[e], domain, DOFs);
		insertElementToDomain(_K, _M, DOFs, Ke, Me, Re, fe, domain);
	}

	Me.resize(0, 0);
	Re.resize(0, 0);
	for (size_t i = 0; i < _mesh->faces().size(); i++) {
		if (_mesh->faces()[i]->inDomain(domain)) {
			processFace(step, matrices, _mesh->faces()[i], Ke, Me, Re, fe, solution);
			fillDOFsIndices(_mesh->faces()[i], domain, DOFs);
			insertElementToDomain(_K, _M, DOFs, Ke, Me, Re, fe, domain);
		}
	}

	for (size_t i = 0; i < _mesh->edges().size(); i++) {
		if (_mesh->edges()[i]->inDomain(domain)) {
			processEdge(step, matrices, _mesh->edges()[i], Ke, Me, Re, fe, solution);
			fillDOFsIndices(_mesh->edges()[i], domain, DOFs);
			insertElementToDomain(_K, _M, DOFs, Ke, Me, Re, fe, domain);
		}
	}

	// TODO: make it direct
	if (matrices & Matrices::K) {
		SparseCSRMatrix<eslocal> csrK = _K;
		_instance->K[domain] = csrK;
		_instance->K[domain].mtype = getMatrixType(step, domain);
	}
	if (matrices & Matrices::M) {
		SparseCSRMatrix<eslocal> csrM = _M;
		_instance->M[domain] = csrM;
		_instance->M[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	}

}

void Physics::assembleMatrix(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe)
{
	updateMatrix(step, matrices, e, Ke, Me, Re, fe, {});
}

void Physics::updateMatrix(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution)
{
	processElement(step, matrices, e, Ke, Me, Re, fe, solution);

	DenseMatrix Ki, Mi, Ri, fi;
	for (size_t i = 0; i < e->filledFaces(); i++) {
		processFace(step, matrices, e->face(i), Ki, Mi, Ri, fi, solution);
		Ke += Ki;
		Me += Mi;
		Re += Ri;
		fe += fi;
	}
	for (size_t i = 0; i < e->filledEdges(); i++) {
		processEdge(step, matrices, e->edge(i), Ki, Mi, Ri, fi, solution);
		Ke += Ki;
		Me += Mi;
		Re += Ri;
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
	DOFs.resize(e->nodes() * pointDOFsOffsets().size());
	for (size_t n = 0, i = 0; n < e->nodes(); n++) {
		for (size_t dof = 0; dof < pointDOFsOffsets().size(); dof++, i++) {
			DOFs[i] = _mesh->nodes()[e->node(n)]->DOFIndex(domain, pointDOFsOffsets()[dof]);
		}
	}
}

void Physics::insertElementToDomain(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, const std::vector<eslocal> &DOFs, const DenseMatrix &Ke, const DenseMatrix &Me, const DenseMatrix &Re, const DenseMatrix &fe, size_t domain)
{
	if (Ke.rows() == DOFs.size() && Ke.columns() == DOFs.size()) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			for (size_t c = 0; c < DOFs.size(); c++) {
				K(DOFs[r], DOFs[c]) = Ke(r, c);
			}
		}
	} else {
		if (Ke.rows() != 0 || Ke.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling stiffness matrix K(" << Ke.rows() << "," << Ke.columns() << ").";
		}
	}

	if (Me.rows() && Me.columns() && DOFs.size() % Me.rows() == 0 && DOFs.size() % Me.columns() == 0) {
		size_t multiplicity = DOFs.size() / Me.rows();
		for (size_t m = 0; m < multiplicity; m++) {
			for (size_t r = 0; r < Me.rows(); r++) {
				for (size_t c = 0; c < Me.columns(); c++) {
					M(DOFs[r * multiplicity + m], DOFs[c * multiplicity + m]) = Me(r, c);
				}
			}
		}
	} else {
		if (Me.rows() != 0 || Me.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling mass matrix M(" << Me.rows() << "," << Me.columns() << ").";
		}
	}


	if (Re.rows() == DOFs.size() && Re.columns() == 1) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			_instance->R[domain][DOFs[r]] += Re(r, 0);
		}
	} else {
		if (Re.rows() != 0 || Re.columns() != 0) {
			std::cout << DOFs;
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling matrix R(" << Re.rows() << "," << Re.columns() << ") with residual forces.";
		}
	}

	if (fe.rows() == DOFs.size() && fe.columns() == 1) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			_instance->f[domain][DOFs[r]] += fe(r, 0);
		}
	} else {
		if (fe.rows() != 0 || fe.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling right-hand side vector f(" << fe.rows() << "," << fe.columns() << ").";
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
		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
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

	// clusterOffset not work because only to the lowest rank receive data
	std::vector<std::vector<std::vector<eslocal> > > skipped(threads, std::vector<std::vector<eslocal> >(_mesh->neighbours().size()));

	std::vector<std::vector<Region*> > allowed(pointDOFs().size());
	if (restriction == SumRestriction::DIRICHLET) {
		allowed = _mesh->getRegionsWithProperties(loadStep, pointDOFs());
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		double tSum = 0;
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {

			if (!_mesh->nodes()[n]->parentElements().size()) {
				// mesh generator can generate dangling nodes -> skip them
				continue;
			}

			if (_mesh->nodes()[n]->clusters().size() > 1) {
				if (_mesh->nodes()[n]->clusters().front() == environment->MPIrank) {
					for (auto c = _mesh->nodes()[n]->clusters().begin() + 1; c != _mesh->nodes()[n]->clusters().end(); ++c) {
						rBufferSize[t][n2i(*c)] += pointDOFs().size();
					}
					incomplete[t].push_back(n);
					incompleteData[t].insert(incompleteData[t].end(), pointDOFs().size(), 0);
					for (auto d = _mesh->nodes()[n]->domains().begin(); d != _mesh->nodes()[n]->domains().end(); ++d) {
						for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
							if (restriction == SumRestriction::NONE || _mesh->commonRegion(allowed[dof], _mesh->nodes()[n]->regions())) {
								eslocal index = _mesh->nodes()[n]->DOFIndex(*d, dof);
								if (index >= 0) {
									incompleteData[t][incompleteData[t].size() - pointDOFs().size() + dof] += data[*d][index];
								}
							}
						}
					}
				} else {
					eslocal cluster = _mesh->nodes()[n]->clusters().front();
					for (auto c = _mesh->nodes()[n]->clusters().begin() + 1; c != _mesh->nodes()[n]->clusters().end(); ++c) {
						if (*c != environment->MPIrank) {
							skipped[t][n2i(*c)].push_back(_mesh->nodes()[n]->clusterOffset(*c));
						}
					}
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
		#pragma omp atomic
		csum += tSum;
	}

	for (size_t n = 0; n < _mesh->neighbours().size(); n++) {
		size_t rSize = rBufferSize[0][n];
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
			skipped[0][n].insert(skipped[0][n].end(), skipped[t][n].begin(), skipped[t][n].end());
			rSize += rBufferSize[t][n];
		}
		rBuffer[n].resize(rSize);
	}

	for (size_t t = 1; t < threads; t++) {
		incomplete[0].insert(incomplete[0].end(), incomplete[t].begin(), incomplete[t].end());
		incompleteData[0].insert(incompleteData[0].end(), incompleteData[t].begin(), incompleteData[t].end());
	}

	#pragma omp parallel for
	for (size_t n = 0; n < _mesh->neighbours().size(); n++) {
		std::sort(skipped[0][n].begin(), skipped[0][n].end());
	}

	if (!Communication::receiveUpperKnownSize(sBuffer[0], rBuffer, _mesh->neighbours())) {
		ESINFO(ERROR) << "problem while exchange sum of squares.";
	}

	auto clusterOffset = [&] (size_t n, eslocal cluster) {
		eslocal offset = _mesh->nodes()[n]->clusterOffset(cluster);
		auto it = std::lower_bound(skipped[0][n2i(cluster)].begin(), skipped[0][n2i(cluster)].end(), offset);
		return offset - (it - skipped[0][n2i(cluster)].begin());
	};

	distribution = Esutils::getDistribution(threads, incomplete[0].size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		double tSum = 0;
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			size_t n = incomplete[0][i];

			for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
				double sum = incompleteData[0][i * pointDOFs().size() + dof];
				for (auto c = _mesh->nodes()[n]->clusters().begin() + 1; c != _mesh->nodes()[n]->clusters().end(); ++c) {
					sum += rBuffer[n2i(*c)][clusterOffset(n, *c) * pointDOFs().size() + dof];
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
		#pragma omp atomic
		csum += tSum;
	}

	MPI_Allreduce(&csum, &gsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return gsum;
}

void Physics::assembleB1(const Step &step, bool withRedundantMultipliers, bool withScaling)
{
	EqualityConstraints::insertDirichletToB1(*_instance, _mesh->regions(), _mesh->nodes(), pointDOFs(), _nodesDOFsOffsets, withRedundantMultipliers);
	EqualityConstraints::insertElementGluingToB1(*_instance, _mesh->neighbours(), _mesh->regions(), _mesh->nodes(), pointDOFs(), _nodesDOFsOffsets, withRedundantMultipliers, withScaling);
}


