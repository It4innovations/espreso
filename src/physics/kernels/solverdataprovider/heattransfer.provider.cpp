
#include "heattransfer.provider.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "math/vector.dense.h"
#include "math/matrix.type.h"
#include "math/matrix.dense.feti.h"
#include "math/matrix.csr.feti.h"

using namespace espreso;

MatrixType HeatTransferSolverDataProvider::General::getMatrixType()
{
	if (	(_configuration.translation_motions.size()) ||
			(_configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR && _configuration.nonlinear_solver.tangent_matrix_correction)
			) {

		return MatrixType::REAL_UNSYMMETRIC;
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void HeatTransferSolverDataProvider::General::dirichletIndices(std::vector<std::pair<esint, esint>> &indices)
{
	for (auto it = _configuration.temperature.begin(); it != _configuration.temperature.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		for (auto i = region->nodes->datatarray().begin(); i != region->nodes->datatarray().end(); ++i) {
			indices.push_back(std::make_pair(*i, 0));
		}
	}
}

void HeatTransferSolverDataProvider::General::dirichletValues(std::vector<double> &values)
{
	size_t offset = 0;
	for (auto it = _configuration.temperature.begin(); it != _configuration.temperature.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		it->second.evaluator->evalSelectedSparse(
				region->nodes->datatarray().size(),
				region->nodes->datatarray().data(),
				Evaluator::Params().coords(3, reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data())),
				values.data() + offset);
		offset += region->nodes->datatarray().size();
	}
}

void HeatTransferSolverDataProvider::General::inequalityIndices(std::vector<std::pair<esint, esint>> &indices)
{

}

void HeatTransferSolverDataProvider::General::inequalityNormals(std::vector<double> &values)
{

}

void HeatTransferSolverDataProvider::General::inequalityGaps(std::vector<double> &values)
{

}

MatrixType HeatTransferSolverDataProvider::FETI::getMatrixType(esint domain)
{
	if (_configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR && _configuration.nonlinear_solver.tangent_matrix_correction) {
		return MatrixType::REAL_UNSYMMETRIC;
	}

	if (_configuration.translation_motions.size()) {
		for (auto it = _configuration.translation_motions.begin(); it != _configuration.translation_motions.end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			for (esint i = info::mesh->elements->eintervalsDistribution[domain]; i < info::mesh->elements->eintervalsDistribution[domain + 1]; i++) {
				if (region->eintervals[i].begin != region->eintervals[i].end) {
					return MatrixType::REAL_UNSYMMETRIC;
				}
			}
		}
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

bool HeatTransferSolverDataProvider::FETI::hasKernel(esint domain)
{
	if (
			_configuration.feti.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
			_configuration.feti.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {
		return true;
	}

	if (_configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT) {
		return false;
	}

	if (_configuration.convection.size()) {
		for (auto it = _configuration.convection.begin(); it != _configuration.convection.end(); ++it) {
			BoundaryRegionStore *region = info::mesh->bregion(it->first);
			if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
				return false;
			}
		}
	}

	if (_configuration.diffuse_radiation.size()) {
		for (auto it = _configuration.diffuse_radiation.begin(); it != _configuration.diffuse_radiation.end(); ++it) {
			BoundaryRegionStore *region = info::mesh->bregion(it->first);
			if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
				return false;
			}
		}
	}

	if (_configuration.bio_heat.size()) {
		for (auto it = _configuration.bio_heat.begin(); it != _configuration.bio_heat.end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			for (esint i = info::mesh->elements->eintervalsDistribution[domain]; i < info::mesh->elements->eintervalsDistribution[domain + 1]; i++) {
				if (region->eintervals[i].begin != region->eintervals[i].end) {
					return false;
				}
			}
		}
	}
	return true;
}

int HeatTransferSolverDataProvider::FETI::initKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	for (esint d = 0; d < K.domains; ++d) {
		if (K[d].type != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
			eslog::error("Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC.\n");
		}
	}

	N1.initDomains(K.domains);
	N2.initDomains(K.domains);
	RegMat.initDomains(K.domains);

	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			N1[d].resize(K[d].nrows, 1);
			N1[d].type = MatrixType::REAL_UNSYMMETRIC;

			RegMat[d].resize(K[d].nrows, K[d].ncols, 1);
			RegMat[d].type = K.type;

			RegMat[d].rows[0] = 1;
			std::fill(RegMat[d].rows + 1, RegMat[d].rows + RegMat[d].nrows + 1, 2);
			RegMat[d].cols[0] = 1;
		}
	}
	return 1;
}

void HeatTransferSolverDataProvider::FETI::fillKernels(MatrixCSRFETI &K, MatrixCSRFETI &M, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			VectorDense diag(K[d].nrows, N1[d].vals);
			K[d].fillDiagonal(&diag);
			RegMat[d].vals[0] = diag.max();

			if (_configuration.feti.conjugate_projector != FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {
				if (ortogonalizeCluster) {
					esint nSum = 0;
					for (esint dd = 0; dd < info::mesh->domains->size; dd++) {
						if (info::mesh->domains->cluster[d] == info::mesh->domains->cluster[dd]) {
							nSum += K[dd].nrows;
						}
					}
					N1[d].fill(1 / std::sqrt(nSum));
				} else {
					N1[d].fill(1 / std::sqrt(K[d].nrows));
				}
			}
		}
	}
}

int HeatTransferSolverDataProvider::Hypre::numfnc()
{
	return 1;
}

void HeatTransferSolverDataProvider::Hypre::initKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N)
{
	// keep empty since the kernel is constant
}

void HeatTransferSolverDataProvider::Hypre::fillKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N)
{
	// keep empty since the kernel is constant
}
