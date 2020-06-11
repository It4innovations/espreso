
#include "structuralmechanics2d.provider.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "basis/utilities/utils.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "config/ecf/input/decomposition.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/envinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/fetidatastore.h"
#include "math/math.h"
#include "math/matrix.type.h"
#include "math/matrix.dense.feti.h"
#include "math/matrix.csr.feti.h"
#include "math/matrix.csr.distributed.h"
#include "math/vector.dense.h"
#include "math/vector.dense.distributed.h"
#include "output/visualization/debug.h"
#include "wrappers/metis/w.metis.h"

#include <cmath>
#include <algorithm>
#include <numeric>

using namespace espreso;

MatrixType StructuralMechanics2DSolverDataProvider::General::getMatrixType()
{
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void StructuralMechanics2DSolverDataProvider::General::dirichletIndices(std::vector<std::pair<esint, esint>> &indices)
{
	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		if (it->second.all.value.size() || it->second.x.value.size()) {
			for (auto i = region->nodes->datatarray().begin(); i != region->nodes->datatarray().end(); ++i) {
				indices.push_back(std::make_pair(*i, 0));
			}
		}
		if (it->second.all.value.size() || it->second.y.value.size()) {
			for (auto i = region->nodes->datatarray().begin(); i != region->nodes->datatarray().end(); ++i) {
				indices.push_back(std::make_pair(*i, 1));
			}
		}
	}
}

void StructuralMechanics2DSolverDataProvider::General::dirichletValues(std::vector<double> &values)
{
	size_t offset = 0;
	double *coors = reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data());
	auto eval = [&] (Evaluator *evaluator, tarray<esint> &nodes) {
		evaluator->evalSelectedSparse(nodes.size(), nodes.data(), Evaluator::Params().coords(3, coors), values.data() + offset);
		offset += nodes.size();
	};

	auto pick = [&] (ECFExpressionOptionalVector &vector, tarray<esint> &nodes) {
		if (vector.all.value.size()) {
			eval(vector.all.evaluator, nodes);
			eval(vector.all.evaluator, nodes);
		} else {
			if (vector.x.value.size()) {
				eval(vector.x.evaluator, nodes);
			}
			if (vector.y.value.size()) {
				eval(vector.y.evaluator, nodes);
			}
		}
	};

	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		pick(it->second, region->nodes->datatarray());
	}
}

void StructuralMechanics2DSolverDataProvider::General::inequalityIndices(std::vector<std::pair<esint, esint>> &indices)
{

}

void StructuralMechanics2DSolverDataProvider::General::inequalityNormals(std::vector<double> &values)
{

}

void StructuralMechanics2DSolverDataProvider::General::inequalityGaps(std::vector<double> &values)
{

}

void addFixPoints(const serializededata<esint, esint>* elements, esint begin, esint end, const serializededata<esint, Element*>* epointers, std::vector<esint> &fixPoints)
{
	esint FIX_POINTS_SIZE = 8;

	auto neighs = [] (std::vector<esint> &neighs, Element::CODE code, int node, const esint* nodes) {
		switch (code) {
		case Element::CODE::HEXA8:
		case Element::CODE::HEXA20:
			if (node < 4) {
				neighs.push_back((nodes[(node + 1) % 4]));
				neighs.push_back((nodes[(node + 3) % 4]));
				neighs.push_back((nodes[node + 4]));
			} else {
				neighs.push_back((nodes[(node + 1) % 4 + 4]));
				neighs.push_back((nodes[(node + 3) % 4 + 4]));
				neighs.push_back((nodes[node - 4]));
			}
			return 3;
		case Element::CODE::TETRA4:
		case Element::CODE::TETRA10:
			neighs.push_back(nodes[(node + 1) % 4]);
			neighs.push_back(nodes[(node + 2) % 4]);
			neighs.push_back(nodes[(node + 3) % 4]);
			return 3;
		case Element::CODE::PRISMA6:
		case Element::CODE::PRISMA15:
			if (node < 3) {
				neighs.push_back(nodes[(node + 1) % 3]);
				neighs.push_back(nodes[(node + 2) % 3]);
				neighs.push_back(nodes[node + 3]);
			} else {
				neighs.push_back(nodes[(node + 1) % 3 + 3]);
				neighs.push_back(nodes[(node + 2) % 3 + 3]);
				neighs.push_back(nodes[node - 3]);
			}
			return 3;

		case Element::CODE::PYRAMID5:
		case Element::CODE::PYRAMID13:
			if (node == 4) {
				neighs.insert(neighs.end(), nodes, nodes + 4);
				return 4;
			} else {
				neighs.push_back(nodes[(node + 1) % 4]);
				neighs.push_back(nodes[(node + 3) % 4]);
				neighs.push_back(nodes[4]);
				return 3;
			}

		case Element::CODE::TRIANGLE3:
		case Element::CODE::TRIANGLE6:
			neighs.push_back(nodes[(node + 1) % 3]);
			neighs.push_back(nodes[(node + 2) % 3]);
			return 2;

		case Element::CODE::SQUARE4:
		case Element::CODE::SQUARE8:
			neighs.push_back(nodes[(node + 1) % 4]);
			neighs.push_back(nodes[(node + 3) % 4]);
			return 2;

		case Element::CODE::LINE2:
		case Element::CODE::LINE3:
			neighs.push_back(nodes[(node + 1) % 2]);
			return 1;
		case Element::CODE::POINT1:
		default:
			return 0;
		}
		return 0;
	};

	std::vector<esint> originnodes, neighsnodes;
	auto element = elements->begin() + begin;
	const auto &epointer = epointers->datatarray();
	for (esint e = 0; e < end - begin; ++e, ++element) {
		for (int n = 0; n < epointer[begin + e]->coarseNodes; ++n) {
			originnodes.insert(
					originnodes.end(),
					neighs(neighsnodes, epointer[begin + e]->code, n, element->data()),
					element->at(n));
		}
	}

	std::vector<esint> permutation(originnodes.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		return originnodes[i] < originnodes[j];
	});

	std::vector<esint> ids, dist, data;
	dist.push_back(0);
	ids.push_back(originnodes[permutation[0]]);
	for (size_t i = 0; i < permutation.size(); i++) {
		if (i && originnodes[permutation[i]] != originnodes[permutation[i - 1]]) {
			utils::sortAndRemoveDuplicates(data, dist.back());
			dist.push_back(data.size());
			ids.push_back(originnodes[permutation[i]]);
		}
		data.push_back(neighsnodes[permutation[i]]);
	}
	utils::sortAndRemoveDuplicates(data, dist.back());
	dist.push_back(data.size());

	for (size_t i = 0; i < data.size(); i++) {
		data[i] = std::lower_bound(ids.begin(), ids.end(), data[i]) - ids.begin();
	}

	std::vector<esint> partition(ids.size());
	METISConfiguration options;
	METIS::call(options, ids.size(), dist.data(), data.data(), 0, NULL, NULL, FIX_POINTS_SIZE, partition.data());

	std::vector<std::vector<esint> > pids(FIX_POINTS_SIZE), pdist(FIX_POINTS_SIZE, { 1 }), pdata(FIX_POINTS_SIZE);
	for (size_t i = 0; i < partition.size(); i++) {
		pids[partition[i]].push_back(i);
	}
	for (size_t i = 0; i < partition.size(); i++) {
		esint p = partition[i];
		for (esint j = dist[i]; j < dist[i + 1]; j++) {
			if (partition[data[j]] == p) {
				size_t index = std::lower_bound(pids[p].begin(), pids[p].end(), data[j]) - pids[p].begin();
				if (pdist[p].size() <= index) {
					pdata[p].push_back(index + 1);
				}
			}
		}
		pdist[p].push_back(pdata[p].size() + 1);
	}

	for (esint p = 0; p < FIX_POINTS_SIZE; p++) {
		if (pids[p].size()) {
			std::vector<float> vals(pdata[p].size(), 1), x(pids[p].size(), 1. / pids[p].size()), y(pids[p].size());
			float last_l = pids[p].size(), l = 1;

			while (fabs((l - last_l) / l) > 1e-6) {
				MATH::upCSRMatVecProduct(pids[p].size(), pids[p].size(), pdist[p].data(), pdata[p].data(), vals.data(), x.data(), y.data());
				last_l = l;
				l = MATH::vecNorm(pids[p].size(), y.data());
				MATH::vecScale(pids[p].size(), 1 / l, y.data());
				x.swap(y);
			}

			fixPoints.push_back(ids[pids[p][MATH::vecNormMaxIndex(pids[p].size(), x.data())]]);
		}
	}
}

StructuralMechanics2DSolverDataProvider::FETI::~FETI()
{
	if (_RegMat) { delete _RegMat; }
}

MatrixType StructuralMechanics2DSolverDataProvider::FETI::getMatrixType(esint domain)
{
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

bool StructuralMechanics2DSolverDataProvider::FETI::hasKernel(esint domain)
{
	if (_configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT) {
		return false;
	}
	return true;
}

int StructuralMechanics2DSolverDataProvider::FETI::initKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	if (info::mesh->FETIData->innerFixPoints.size() == 0) {

	}

	for (esint d = 0; d < K.domains; ++d) {
		if (K[d].type != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
			eslog::error("Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC.\n");
		}
	}

	size_t clusters = *std::max_element(info::mesh->elements->clusters.begin(), info::mesh->elements->clusters.end()) + 1;

	_cCenter = _cNorm = std::vector<Point>(clusters, Point(0, 0, 0));
	_cNp = std::vector<size_t>(clusters, 0);

	_dCenter = _dNorm = std::vector<Point>(info::mesh->elements->ndomains, Point(0, 0, 0));

	std::vector<double> cbuffer1(info::mesh->elements->ndomains, 0);
	dnodes.resize(info::mesh->elements->ndomains);

	auto dmap = info::mesh->nodes->domains->cbegin();
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				_dCenter[domain] += info::mesh->nodes->coordinates->datatarray()[i];
				dnodes[domain].push_back(i);
			}
		}
	}

	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		_cCenter[info::mesh->elements->clusters[d]] += _dCenter[d];
		_dCenter[d] = _dCenter[d] / dnodes[d].size();
		_cNp[info::mesh->elements->clusters[d]] += dnodes[d].size();
	}
	for (size_t c = 0; c < clusters; c++) {
		_cCenter[c] /= _cNp[c];
	}

	// Compute norm of column 4 (norm.x)
	dmap = info::mesh->nodes->domains->cbegin();
	std::vector<double> pnorm(info::mesh->elements->ndomains), pcnorm(info::mesh->elements->ndomains);
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				Point dp = info::mesh->nodes->coordinates->datatarray()[i] - _dCenter[domain];
				pnorm[domain] += dp.x * dp.x + dp.y * dp.y;
				Point cp = info::mesh->nodes->coordinates->datatarray()[i] - _cCenter[info::mesh->elements->clusters[domain]];
				pcnorm[domain] += cp.x * cp.x + cp.y * cp.y;
			}
		}
	}
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		_dNorm[d].x = std::sqrt(pnorm[d]);
		cbuffer1[d] += pcnorm[d];
	}
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		_cNorm[info::mesh->elements->clusters[d]].x += cbuffer1[d];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].x = std::sqrt(_cNorm[c].x);
	}

	N1.initDomains(K.domains);
	N2.initDomains(K.domains);
	RegMat.initDomains(K.domains);
	_RegMat = new MatrixCSRFETI();
	_RegMat->initDomains(info::mesh->elements->ndomains);

	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			std::vector<esint> fixPoints;
		//	if (_BEMDomain[domain]) {
		//		fixPoints = std::vector<esint>(
		//				info::mesh->FETIData->surfaceFixPoints.begin() + info::mesh->FETIData->sFixPointsDistribution[domain],
		//				info::mesh->FETIData->surfaceFixPoints.begin() + info::mesh->FETIData->sFixPointsDistribution[domain + 1]);
		//	} else {
				fixPoints = std::vector<esint>(
						info::mesh->FETIData->innerFixPoints.begin() + info::mesh->FETIData->iFixPointsDistribution[d],
						info::mesh->FETIData->innerFixPoints.begin() + info::mesh->FETIData->iFixPointsDistribution[d + 1]);
		//	}

			MatrixCSR Nt(3, K[d].ncols, 4 * fixPoints.size());

			Nt.rows[0] = 1;
			Nt.rows[1] = Nt.rows[0] + fixPoints.size();
			Nt.rows[2] = Nt.rows[1] + fixPoints.size();
			Nt.rows[3] = Nt.rows[2] + 2 * fixPoints.size();

			auto n2DOF = [&] (esint node) {
				return std::lower_bound(dnodes[d].begin(), dnodes[d].end(), node) - dnodes[d].begin();
			};

			esint cindex = 0;
			for (size_t c = 0; c < 2; c++) {
				for (size_t i = 0; i < fixPoints.size(); i++, cindex++) {
					Nt.cols[cindex] = 2 * n2DOF(fixPoints[i]) + c + 1;
					Nt.vals[cindex] = 1;
				}
			}

			for (size_t i = 0; i < fixPoints.size(); i++, cindex += 2) {
				const Point &p = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
				Nt.cols[cindex]     = 2 * n2DOF(fixPoints[i]) + 0 + 1;
				Nt.cols[cindex + 1] = 2 * n2DOF(fixPoints[i]) + 1 + 1;
				Nt.vals[cindex]     = -p.y;
				Nt.vals[cindex + 1] =  p.x;
			}

			// N * (tran(N) * N)^-1 * tran(N)
			//
			// AX = B => X = A^-1B
			//
			// if A = tran(N) * N;
			// then X = (tran(N) * N)^-1 * N
			// RegMat = tran(N) * X

			MatrixCSR N;
			Nt.transposeTo(&N);
			MatrixCSR A;

			A.multiply(Nt, N);
			A.removeLower(MatrixType::REAL_SYMMETRIC_INDEFINITE);
			MatrixDense _X(N.nrows, N.ncols), B = N;
			A.solve(B, _X);
			MatrixCSR X = _X;
			X.transpose();
			_RegMat->at(d)->multiply(N, X);
			_RegMat->at(d)->removeLower(MatrixType::REAL_SYMMETRIC_INDEFINITE);

			RegMat[d].shallowCopyStructure(_RegMat->at(d));

			// Gram-schmidt

			N1[d].resize(K[d].nrows, 3);
		}
	}
	return 6;
}

void StructuralMechanics2DSolverDataProvider::FETI::fillKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			VectorDense diag(K[d].nrows, N1[d].vals);
			K[d].fillDiagonal(&diag);

			RegMat[d].fillData(_RegMat->at(d));
			RegMat[d].scale(diag.max());

			Point center = _dCenter[d], norm = _dNorm[d];
			size_t np = dnodes[d].size();

			if (ortogonalizeCluster) {
				size_t cluster = info::mesh->elements->clusters[d];
				center = _cCenter[cluster], norm = _cNorm[cluster];
				np = _cNp[cluster];
			} else {
				center = _dCenter[d], norm = _dNorm[d];
				np = dnodes[d].size();
			}

			double v = 1 / std::sqrt(np);
			for (size_t i = 0; i < dnodes[d].size(); ++i) {
				Point p = info::mesh->nodes->coordinates->datatarray()[dnodes[d][i]] - center;

				N1[d][2 * i + 0][0] = v;
				N1[d][2 * i + 1][0] = 0;

				N1[d][2 * i + 0][1] = 0;
				N1[d][2 * i + 1][1] = v;

				N1[d][2 * i + 0][2] = -p.y / norm.x;
				N1[d][2 * i + 1][2] =  p.x / norm.x;
			}
		}
	}
}

int StructuralMechanics2DSolverDataProvider::Hypre::numfnc()
{
	return 2;
}

void StructuralMechanics2DSolverDataProvider::Hypre::initKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N)
{
	N.initVectors(3);
	N.resize(K.nrows, K.nhalo, K.nneighbors);
}

void StructuralMechanics2DSolverDataProvider::Hypre::fillKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N)
{
	for (esint n = 0; n < info::mesh->nodes->size; n++) {
		Point p = info::mesh->nodes->coordinates->datatarray()[n];
		N[0][2 * n + 0] = 1;
		N[0][2 * n + 1] = 0;

		N[1][2 * n + 0] = 0;
		N[1][2 * n + 1] = 1;

		N[2][2 * n + 0] =  p.y;
		N[2][2 * n + 1] = -p.x;
	}
}
