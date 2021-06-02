
#include "provider.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "config/ecf/input/decomposition.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "math/math.h"
#include "math/matrix.ijv.feti.h"
#include "math/matrix.csr.feti.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/fetidatastore.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "output/visualization/debug.h"
#include "wrappers/metis/w.metis.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void SolverDataProvider::FETI::computeCornerNodes()
{
	if (info::mesh->nodes->domains == NULL) {
		mesh::computeNodeDomainDistribution(info::mesh->elements, info::mesh->nodes, info::mesh->domains, info::mesh->neighborsWithMe);
	}

	std::vector<esint> uniq;
	std::vector<std::vector<esint> > nodes;

	esint i = 0;
	for (auto dmap = info::mesh->nodes->domains->cbegin(); dmap != info::mesh->nodes->domains->cend(); ++dmap, ++i) {
		if (dmap->size() > 1) {
			size_t index = 0;
			for (auto u = uniq.begin(); u != uniq.end(); ++u, ++index) {
				if (*u == (esint)dmap->size() && memcmp(&*(u + 1), dmap->data(), dmap->size() * sizeof(esint)) == 0) {
					break;
				}
				u += *u;
			}
			if (index == nodes.size()) {
				uniq.push_back(dmap->size());
				uniq.insert(uniq.end(), dmap->begin(), dmap->end());
				nodes.push_back({ i });
			} else {
				nodes[index].push_back(i);
			}
		}
	}

	for (size_t n = 0; n < nodes.size(); ++n) {
		esint inc = nodes[n].size() / 3;
		inc = inc ? inc : 1;
		for (size_t nn = 0; nn < nodes[n].size(); nn += inc) {
			info::mesh->FETIData->corners.push_back(nodes[n][nn]);
		}
	}

	std::sort(info::mesh->FETIData->corners.begin(), info::mesh->FETIData->corners.end());

	eslog::checkpointln("MESH: CORNER NODES COMPUTED");
}

static void addFixPoints(const serializededata<esint, esint>* elements, esint begin, esint end, const serializededata<esint, Element*>* epointers, std::vector<esint> &fixPoints)
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

void SolverDataProvider::FETI::computeFixPoints()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > fixPoints(threads), fixPointsDist(threads);

	fixPointsDist.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> dist, data, partition;
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			size_t size = fixPoints[t].size();
			addFixPoints(info::mesh->elements->nodes, info::mesh->domains->elements[d], info::mesh->domains->elements[d + 1], info::mesh->elements->epointers, fixPoints[t]);
			utils::sortAndRemoveDuplicates(fixPoints[t], size);
			fixPointsDist[t].push_back(fixPoints[t].size());
		}
	}

	utils::threadDistributionToFullDistribution(fixPointsDist);
	info::mesh->FETIData->iFixPointsDistribution.clear();
	for (size_t t = 0; t < threads; t++) {
		info::mesh->FETIData->iFixPointsDistribution.insert(info::mesh->FETIData->iFixPointsDistribution.end(), fixPointsDist[t].begin(), fixPointsDist[t].end());
		info::mesh->FETIData->innerFixPoints.insert(info::mesh->FETIData->innerFixPoints.end(), fixPoints[t].begin(), fixPoints[t].end());
	}

	DebugOutput::surfaceFixPoints();
}

void SolverDataProvider::FETI::computeFixPointsOnSurface()
{
	if (info::mesh->domainsSurface == NULL || info::mesh->domainsSurface->enodes == NULL) {
		mesh::computeDomainsSurface(info::mesh->nodes, info::mesh->elements, info::mesh->domains, info::mesh->domainsSurface, info::mesh->neighbors);
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > fixPoints(threads), fixPointsDist(threads);

	fixPointsDist.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> dist, data, partition;
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			size_t size = fixPoints[t].size();
			addFixPoints(info::mesh->domainsSurface->enodes, info::mesh->domainsSurface->edistribution[d], info::mesh->domainsSurface->edistribution[d + 1], info::mesh->domainsSurface->epointers, fixPoints[t]);
			utils::sortAndRemoveDuplicates(fixPoints[t], size);
			fixPointsDist[t].push_back(fixPoints[t].size());
		}
	}

	utils::threadDistributionToFullDistribution(fixPointsDist);
	info::mesh->FETIData->sFixPointsDistribution.clear();
	for (size_t t = 0; t < threads; t++) {
		info::mesh->FETIData->sFixPointsDistribution.insert(info::mesh->FETIData->sFixPointsDistribution.end(), fixPointsDist[t].begin(), fixPointsDist[t].end());
		info::mesh->FETIData->surfaceFixPoints.insert(info::mesh->FETIData->surfaceFixPoints.end(), fixPoints[t].begin(), fixPoints[t].end());
	}

	DebugOutput::surfaceFixPoints();
}

void SolverDataProvider::FETI::_buildB0FromCorners(MatrixCSRFETI &K, MatrixIJVFETI &B0, int DOFs)
{
	if (info::mesh->FETIData->corners.size() == 0) {

	}

	esint lambda = 0;
	std::vector<std::vector<esint> > ROWS(info::mesh->domains->size), COLS(info::mesh->domains->size);
	std::vector<std::vector<double> > VALS(info::mesh->domains->size);

	auto dmap = K.dmap->begin();
	for (size_t n = 0, prev = 0; n < info::mesh->FETIData->corners.size(); prev = n++) {
		dmap += DOFs * (info::mesh->FETIData->corners[n] - info::mesh->FETIData->corners[prev]);
		for (int dof = 0; dof < DOFs; ++dof) {
			for (auto di1 = (dmap + dof)->begin(), di2 = di1 + 1; di2 != (dmap + dof)->end(); ++di1, ++di2) {
				if (K.ismy(di1->domain) && K.ismy(di2->domain)) {
					ROWS[di1->domain - K.doffset].push_back(lambda + 1);
					COLS[di1->domain - K.doffset].push_back(di1->index + 1);
					VALS[di1->domain - K.doffset].push_back(1);

					ROWS[di2->domain - K.doffset].push_back(lambda + 1);
					COLS[di2->domain - K.doffset].push_back(di2->index + 1);
					VALS[di2->domain - K.doffset].push_back(-1);
					++lambda;
				}
			}
		}
	}

	B0.initDomains(K.domains);
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->domains->size; d++) {
		B0[d].type = MatrixType::REAL_UNSYMMETRIC;
		B0[d].resize(lambda, K[d].ncols, VALS[d].size());
		B0[d].fillPattern(ROWS[d].size(), ROWS[d].data(), COLS[d].data());
		B0[d].fillValues(VALS[d].size(), VALS[d].data());
	}
}
