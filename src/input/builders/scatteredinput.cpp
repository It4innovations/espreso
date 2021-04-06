
#include "scatteredinput.h"
#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/communication.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

ScatteredInput::ScatteredInput(MeshBuilder &meshData)
: Input(meshData), _sfc(info::mesh->dimension, SFCDEPTH, _meshData.coordinates.size(), _meshData.coordinates.data())
{
	profiler::syncstart("scattered_input");

	if (info::mpi::size == 1) {
		eslog::globalerror("ESPRESO internal error: use the sequential input for building mesh on 1 MPI process.\n");
	}

	eslog::startln("BUILDER: BUILD SCATTERED MESH", "BUILDER");

	balance();
	eslog::checkpointln("BUILDER: DATA BALANCED");

	assignRegions(_meshData.eregions, _meshData.eIDs, _eDistribution, _eregsize, _eregions);
	assignRegions(_meshData.nregions, _meshData.nIDs, _nDistribution, _nregsize, _nregions);
	eslog::checkpointln("BUILDER: REGION ASSIGNED");

//	reindexRegions();

	assignNBuckets();
	profiler::synccheckpoint("nbuckets_compputed");
	eslog::checkpointln("BUILDER: NODES BUCKETS COMPUTED");

	assignEBuckets();
	profiler::synccheckpoint("ebuckets_compputed");
	eslog::checkpointln("BUILDER: ELEMENTS BUCKETS ASSIGNED");

	clusterize();
	profiler::synccheckpoint("cluterization");
	eslog::checkpointln("BUILDER: ELEMENTS CLUSTERED");

	computeSFCNeighbors();
	profiler::synccheckpoint("approximate_neighbors");
	eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");

	if (_meshData.removeDuplicates) {
		mergeDuplicatedNodes();
		profiler::synccheckpoint("merge_duplicate_nodes");
		eslog::checkpointln("BUILDER: DUPLICATED NODES MERGED");
	}

	sortElements();
	profiler::synccheckpoint("sort_nodes");
	eslog::checkpointln("BUILDER: ELEMENTS SORTED");

	linkup();
	fillNeighbors();
	profiler::synccheckpoint("linkup");
	eslog::checkpointln("BUILDER: LINKED UP");

	sortNodes();
	profiler::synccheckpoint("final_sort");
	eslog::checkpointln("BUILDER: NODES SORTED");

	reindexElementNodes();
	profiler::synccheckpoint("reindex_elements");
	eslog::checkpointln("BUILDER: ELEMENTS NODES REINDEXED");

	if (_meshData.removeDuplicates) {
		removeDuplicateElements();
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS REMOVED");
		profiler::synccheckpoint("merge_duplicate_elements");
	}

	fillNodes();
	profiler::synccheckpoint("fill_nodes");
	eslog::checkpointln("BUILDER: NODES FILLED");

	fillElements();
	profiler::synccheckpoint("fill_elements");
	eslog::checkpointln("BUILDER: ELEMENTS SORTED");

	if (info::mesh->nodes->elements == NULL) {
		mesh::linkNodesAndElements();
		profiler::synccheckpoint("link_nodes_and_elements");
	}

	exchangeBoundary();
	profiler::synccheckpoint("exchange_boundary");
	eslog::checkpointln("BUILDER: BOUNDARY EXCHANGED");

	fillRegions(_meshData.eregions, _eregsize, _eregions);
	fillRegions(_meshData.nregions, _nregsize, _nregions);
	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	profiler::synccheckpoint("fill_regions");
	eslog::checkpointln("BUILDER: REGIONS FILLED");

	reindexBoundaryNodes();
	profiler::synccheckpoint("reindex_boundary_nodes");
	eslog::endln("BUILDER: BOUNDARY NODES REINDEXED");

	profiler::syncend("scattered_input");
//	polish();
}

void ScatteredInput::assignNBuckets()
{
	profiler::syncstart("assign_nbuckets");
	profiler::syncparam("size", _meshData.nIDs.size());
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<size_t> cdistribution = tarray<size_t>::distribute(threads, _meshData.nIDs.size());

	_nBuckets.resize(_meshData.coordinates.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = cdistribution[t]; n < cdistribution[t + 1]; ++n) {
			_nBuckets[n] = _sfc.getBucket(_meshData.coordinates[n]);
		}
	}
	profiler::syncend("assign_nbuckets");
}

void ScatteredInput::assignEBuckets()
{ // we needs to ask a neighbor process to get bucket of (arbitrary) node -- now the closest process is chosen
	profiler::syncstart("assign_ebuckets");
	profiler::syncparam("size", _meshData.eIDs.size());
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, _meshData.esize.size());

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	std::vector<esint> closest(_meshData.esize.size());
	_eBuckets.resize(_meshData.esize.size());

	esint nbegin = _nDistribution[info::mpi::rank];
	esint nend = _nDistribution[info::mpi::rank + 1];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = edistribution[t], offset = _meshData._edist[edistribution[t]]; e < edistribution[t + 1]; offset += _meshData.esize[e++]) {
			closest[e] = _meshData.enodes[offset];
			for (esint n = 1; n < _meshData.esize[e]; n++) {
				if (nbegin <= _meshData.enodes[offset + n] && _meshData.enodes[offset + n] < nend) {
					if (closest[e] > _meshData.enodes[offset + n] || closest[e] < nbegin) {
						closest[e] = _meshData.enodes[offset + n];
					}
				} else {
					if (std::abs(closest[e] - nbegin) > std::abs(_meshData.enodes[offset + n] - nbegin)) {
						closest[e] = _meshData.enodes[offset + n];
					}
				}
			}
		}
	}

	profiler::synccheckpoint("compute_closest");

	std::vector<esint> permutation(_meshData.esize.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return closest[i] < closest[j]; });

	profiler::synccheckpoint("sort");

	std::vector<esint> sNodes, rNodes;
	std::vector<esint> sBuckets, rBuckets;
	std::vector<int> targets, sources;

	sNodes.reserve(permutation.size() + 2 * info::mpi::size);

	size_t prevsize;
	auto begin = permutation.begin();
	for (int r = 0; r < info::mpi::size; r++) {
		prevsize = sNodes.size();
		sNodes.push_back(0);
		sNodes.push_back(r);
		sNodes.push_back(info::mpi::rank);

		auto n = begin;
		for ( ; n != permutation.end() && closest[*n] < _nDistribution[r + 1]; ++n) {
			sNodes.push_back(closest[*n]);
		}
		sNodes[prevsize] = 3 + n - begin;
		begin = n;
	}
	profiler::synccheckpoint("unknown_nodes");

	if (!Communication::allToAllWithDataSizeAndTarget(sNodes, rNodes)) {
		eslog::error("ESPRESO internal error: ask neighbors for nodes buckets.\n");
	}
	profiler::synccheckpoint("exchange");

	std::vector<esint> boundaries(info::mpi::size);
	size_t offset = 0;
	for (int r = 1; r < info::mpi::size; r++, offset += rNodes[offset]) {
		boundaries[r] = rNodes[offset] + boundaries[r - 1];
	}
	std::sort(boundaries.begin(), boundaries.end(), [&] (esint i, esint j) {
		return rNodes[i + 2] < rNodes[j + 2];
	});

	sBuckets.reserve(rNodes.size());

	for (int r = 0; r < info::mpi::size; r++) {
		offset = boundaries[r];
		esint size = rNodes[offset++] - 3;
		offset++; //skip rank
		offset++; //skip target

		sBuckets.push_back(3 + size);
		sBuckets.push_back(r);
		sBuckets.push_back(info::mpi::rank);
		auto it = _meshData.nIDs.begin();
		for (esint n = 0; n < size; ++n, ++offset) {
			while (*it < rNodes[offset]) { ++it; }
			sBuckets.push_back(_nBuckets[it - _meshData.nIDs.begin()]);
		}
	}
	profiler::synccheckpoint("unknown_nodes_buckets");

	if (!Communication::allToAllWithDataSizeAndTarget(sBuckets, rBuckets)) {
		eslog::error("ESPRESO internal error: return nodes buckets.\n");
	}
	profiler::synccheckpoint("exchange");

	boundaries[0] = offset = 0;
	for (int r = 1; r < info::mpi::size; r++, offset += rBuckets[offset]) {
		boundaries[r] = rBuckets[offset] + boundaries[r - 1];
	}
	std::sort(boundaries.begin(), boundaries.end(), [&] (esint i, esint j) {
		return rBuckets[i + 2] < rBuckets[j + 2];
	});

	size_t e = 0;
	for (int r = 0; r < info::mpi::size; r++) {
		offset = boundaries[r];
		esint size = rBuckets[offset++] - 3;
		offset++; //skip rank
		offset++; //skip target

		for (esint n = 0; n < size; ++n, ++offset) {
			_eBuckets[permutation[e++]] = rBuckets[offset];
		}
	}
	profiler::synccheckpoint("set buckets");
	profiler::syncend("assign_ebuckets");
}

void ScatteredInput::clusterize()
{
	profiler::syncstart("clusterize");
	profiler::syncparam("nbuckets", _nBuckets.size());
	profiler::syncparam("ebuckets", _eBuckets.size());
	size_t threads = info::env::OMP_NUM_THREADS;

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	std::vector<esint> npermutation(_nBuckets.size()), epermutation(_eBuckets.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (esint i, esint j) { return _nBuckets[i] < _nBuckets[j]; });
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return _eBuckets[i] < _eBuckets[j]; });

//	if (!Communication::computeSFCBalancedBorders(_sfc, _eBuckets, epermutation, _bucketsBorders)) {
//		eslog::error("MESIO internal error: cannot balance SFC.\n");
//	}
	profiler::synccheckpoint("sort");

	if (!Communication::computeSplitters(_eBuckets, epermutation, _bucketsBorders)) {
		eslog::error("MESIO internal error: cannot balance SFC.\n");
	}

	profiler::synccheckpoint("compute_splitters");
	_bucketsBorders.back() = _sfc.buckets(_sfc.depth());

	_nregsize = _meshData.nregions.size() / (8 * sizeof(esint)) + 1;
	_eregsize = _meshData.eregions.size() / (8 * sizeof(esint)) + 1;

	_nregions.resize(_nregsize * _meshData.nIDs.size());
	_eregions.resize(_eregsize * _meshData.eIDs.size());

	size_t r = 0;
	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));

		std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, nregion->second.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = rdistribution[t]; i < rdistribution[t + 1]; ++i) {
				_nregions[_nregsize * nregion->second[i] + byte] |= bit;
			}
		}
	}
	r = 0;
	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));

		std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, eregion->second.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = rdistribution[t]; i < rdistribution[t + 1]; ++i) {
				_eregions[_eregsize * eregion->second[i] + byte] |= bit;
			}
		}
	}

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(
			5 * info::mpi::size +
			// esize, eID, etype, body, material, regions
			(5 + _eregsize) * _meshData.esize.size() +
			_meshData.enodes.size() +
			(1 + _nregsize) * _meshData.nIDs.size() +
			_meshData.coordinates.size() * sizeof(Point) / sizeof(esint));

	size_t prevsize;
	auto nbegin = npermutation.begin();
	auto ebegin = epermutation.begin();
	for (int r = 0; r < info::mpi::size; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of elements
		sBuffer.push_back(0); // number of elements nodes
		sBuffer.push_back(0); // number of coordinates

		auto e = ebegin;
		for ( ; e != epermutation.end() && _eBuckets[*e] < _bucketsBorders[r + 1]; ++e) {
			sBuffer.push_back(_meshData.esize[*e]);
			sBuffer.push_back(_meshData.eIDs[*e]);
			sBuffer.push_back(_meshData.etype[*e]);
			sBuffer.push_back(_meshData.body[*e]);
			sBuffer.push_back(_meshData.material[*e]);
			sBuffer.insert(sBuffer.end(), _eregions.begin() + *e * _eregsize, _eregions.begin() + (*e + 1) * _eregsize);
			sBuffer.insert(sBuffer.end(), _meshData.enodes.begin() + _meshData._edist[*e], _meshData.enodes.begin() + _meshData._edist[*e + 1]);
			sBuffer[prevsize + 3] += _meshData._edist[*e + 1] - _meshData._edist[*e];
		}
		sBuffer[prevsize + 2] = e - ebegin;
		ebegin = e;

		auto n = nbegin;
		for ( ; n != npermutation.end() && _nBuckets[*n] < _bucketsBorders[r + 1]; ++n) {
			sBuffer.push_back(_meshData.nIDs[*n]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n + 1));
			sBuffer.insert(sBuffer.end(), _nregions.begin() + *n * _nregsize, _nregions.begin() + (*n + 1) * _nregsize);
		}
		sBuffer[prevsize + 4] = n - nbegin;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("ESPRESO internal error: distribute elements according to SFC.\n");
	}
	profiler::synccheckpoint("exchange");

	_meshData.esize.clear();
	_meshData.eIDs.clear();
	_meshData.etype.clear();
	_meshData.body.clear();
	_meshData.material.clear();
	_meshData.enodes.clear();
	_meshData._edist.clear();
	_eregions.clear();

	_meshData.nIDs.swap(_nIDs); // keep for later usage in linkup phase
	_meshData.coordinates.clear();
	_nregions.clear();

	size_t offset = 0;
	Point point;
	for (int r = 0; r < info::mpi::size; r++) {
		++offset;
		size_t esize = rBuffer[++offset];
		++offset;
		size_t csize = rBuffer[++offset]; // coordinates
		++offset;

		for (size_t e = 0; e < esize; ++e) {
			_meshData.esize.push_back(rBuffer[offset++]);
			_meshData.eIDs.push_back(rBuffer[offset++]);
			_meshData.etype.push_back(rBuffer[offset++]);
			_meshData.body.push_back(rBuffer[offset++]);
			_meshData.material.push_back(rBuffer[offset++]);
			_eregions.insert(_eregions.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _eregsize);
			offset += _eregsize;
			_meshData.enodes.insert(_meshData.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _meshData.esize.back());
			offset += _meshData.esize.back();
		}
		for (size_t c = 0; c < csize; ++c) {
			_meshData.nIDs.push_back(rBuffer[offset]);
			++offset;
			memcpy(reinterpret_cast<void*>(&point), rBuffer.data() + offset, sizeof(Point));
			_meshData.coordinates.push_back(point);
			offset += sizeof(Point) / sizeof(esint);
			_nregions.insert(_nregions.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _nregsize);
			offset += _nregsize;
		}
	}
	profiler::synccheckpoint("rbuffer");

//	if (!_meshData.eIDs.size()) {
//		eslog::error("ESPRESO internal error: a process without elements -- re-run with smaller number of MPI.\n");
//	}
	auto back = _eDistribution.back();
	if (_meshData.eIDs.size()) {
		back = _meshData.eIDs.back();
	}
	Communication::allGather(&back, _eDistribution.data() + 1, sizeof(back), MPI_BYTE);
	for (size_t i = 1; i < _eDistribution.size(); i++) {
		++_eDistribution[i];
	}
	profiler::synccheckpoint("edistribution");
	profiler::syncend("clusterize");
}

void ScatteredInput::computeSFCNeighbors()
{
	profiler::syncstart("sfc_neighbors");
	std::vector<std::pair<size_t, size_t> > neighbors;

	size_t index = _bucketsBorders[info::mpi::rank];
	size_t last = _bucketsBorders[info::mpi::rank + 1];
	while (index < last) {
		size_t depth = _sfc.depth(), bsize = 1;
		while (depth > 1 && index % (bsize * _sfc.bucketSize()) == 0 && index + (bsize * _sfc.bucketSize()) < last) {
			--depth;
			bsize *= _sfc.bucketSize();
		}
		_sfc.addSFCNeighbors(depth, index / bsize, _bucketsBorders, neighbors);
		index += bsize;
	}

	profiler::synccheckpoint("add_sfc_neighbors");
//	_sfc.iterateBuckets(_bucketsBorders[info::mpi::rank], _bucketsBorders[info::mpi::rank + 1], [&] (size_t depth, size_t index) {
//		_sfc.addSFCNeighbors(depth, index, neighbors);
//	});

	utils::sortAndRemoveDuplicates(neighbors);

	for (size_t i = 0; i < neighbors.size(); i++) {
		size_t bstep = _sfc.buckets(_sfc.depth()) / _sfc.buckets(neighbors[i].first);
		neighbors[i].first = neighbors[i].second * bstep;
		neighbors[i].second = neighbors[i].second * bstep + bstep;
	}

	std::sort(neighbors.begin(), neighbors.end());

	// we assume that where is not any interval accross processes (assured by the above addNeighbors alg.)
	auto rank = _bucketsBorders.begin();
	for (size_t i = 0; i < neighbors.size(); i++) {
		while (*rank <= (esint)neighbors[i].first) { ++rank; }
		int r = rank - _bucketsBorders.begin() - 1;
		if (r != info::mpi::rank && (_sfcNeighbors.empty() || _sfcNeighbors.back() != r)) {
			_sfcNeighbors.push_back(r);
		}
	}
	profiler::synccheckpoint("sort_and_unique");
	profiler::syncend("sfc_neighbors");
}

void ScatteredInput::mergeDuplicatedNodes()
{
	profiler::syncstart("merge_duplicated_nodes");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<size_t> cdistribution = tarray<size_t>::distribute(threads, _meshData.nIDs.size());

	double eps = info::ecf->input.duplication_tolerance;
	eps += info::ecf->input.duplication_tolerance * 1e-6;

	std::vector<std::vector<std::pair<int, size_t> > > toneighs(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<int, size_t> > tneighs;
		std::vector<int> neighs;
		neighs.reserve(27);
		for (size_t n = cdistribution[t]; n < cdistribution[t + 1]; ++n) {
			int hit[6] = { 0, 0, 0, 0, 0, 0};
			for (int d = 0; d < info::mesh->dimension; ++d) {
				size_t origin = _sfc.getBucket(_meshData.coordinates[n][d], d);
				if (_sfc.getBucket(_meshData.coordinates[n][d] - eps, d) < origin) {
					hit[2 * d + 0] = -1;
				}
				if (_sfc.getBucket(_meshData.coordinates[n][d] + eps, d) > origin) {
					hit[2 * d + 1] = +1;
				}
			}

			neighs.clear();
			for (int x = hit[0]; x <= hit[1]; ++x) {
				for (int y = hit[2]; y <= hit[3]; ++y) {
					for (int z = hit[4]; z <= hit[5]; ++z) {
						if (x || y || z) {
							size_t b = _sfc.getBucket(_meshData.coordinates[n] + Point(x * eps, y * eps, z * eps));
							int nn = std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), b + 1) - _bucketsBorders.begin() - 1;
							if (nn != info::mpi::rank) {
								if (neighs.size() == 0 || neighs.back() != nn) {
									neighs.push_back(nn);
								}
							}
						}
					}
				}
			}
			utils::sortAndRemoveDuplicates(neighs);
			for (auto nn = neighs.begin(); nn != neighs.end(); ++nn) {
				tneighs.push_back(std::make_pair(*nn, n));
			}
		}
		toneighs[t].swap(tneighs);
	}

	profiler::synccheckpoint("compute_border_nodes");

	std::vector<esint> offsets, ids, regions;
	std::vector<Point> coordinates;
	std::vector<std::vector<esint> > sBuffer(_sfcNeighbors.size()), rBuffer(_sfcNeighbors.size());
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = 0; i < toneighs[t].size(); ++i) {
			size_t ni = std::lower_bound(_sfcNeighbors.begin(), _sfcNeighbors.end(), toneighs[t][i].first) - _sfcNeighbors.begin();
			size_t node = toneighs[t][i].second;

			offsets.push_back(node);
			ids.push_back(_meshData.nIDs[node]);
			coordinates.push_back(_meshData.coordinates[node]);
			for (size_t r = 0; r < _nregsize; ++r) {
				regions.push_back(_nregions[_nregsize * node + r]);
			}

			sBuffer[ni].push_back(_meshData.nIDs[node]);
			esint *begin = reinterpret_cast<esint*>(_meshData.coordinates.data() + node);
			esint *end = reinterpret_cast<esint*>(_meshData.coordinates.data() + node + 1);
			sBuffer[ni].insert(sBuffer[ni].end(), begin, end);
			for (size_t r = 0; r < _nregsize; ++r) {
				sBuffer[ni].push_back(_nregions[_nregsize * node + r]);
			}
		}
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: cannot exchange potentially duplicated nodes.\n");
	}
	profiler::synccheckpoint("exchange");

	std::vector<esint> noffset;
	noffset.push_back(ids.size());
	for (size_t n = 0; n < _sfcNeighbors.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ) {
			ids.push_back(rBuffer[n][i++]);
			Point p;
			memcpy(&p, rBuffer[n].data() + i, sizeof(Point));
			coordinates.push_back(p);
			i += sizeof(p) / sizeof(esint);
			for (size_t r = 0; r < _nregsize; ++r) {
				regions.push_back(rBuffer[n][i++]);
			}
		}
		noffset.push_back(ids.size());
	}
	profiler::synccheckpoint("rbuffer");

	std::vector<size_t> toinsert, toerase;
	searchDuplicateNodes(coordinates, ids, [&] (esint id, esint target) {
		if (id < noffset.front() && noffset.front() <= target) { // I have duplicate
			toerase.push_back(offsets[id]);
		}
		if (target < noffset.front() && noffset.front() <= id) { // rank n has duplicate and I have target
			toinsert.push_back(id);
		}
	});
	profiler::synccheckpoint("search");

	if (toerase.size()) {
		utils::sortAndRemoveDuplicates(toerase);
		for (size_t n = toerase.front(), last = toerase.front(), e = 0; n < _meshData.nIDs.size(); ++n) {
			if (toerase.size() <= e || n != toerase[e]) {
				_meshData.nIDs[last] = _meshData.nIDs[n];
				_meshData.coordinates[last] = _meshData.coordinates[n];
				for (size_t r = 0; r < _nregsize; ++r) {
					_nregions[last * _nregsize + r] = _nregions[n * _nregsize + r];
				}
				++last;
			} else {
				++e;
			}
		}
		_nregions.resize(_nregsize * (_meshData.nIDs.size() - toerase.size()));
		_meshData.coordinates.resize(_meshData.nIDs.size() - toerase.size());
		_meshData.nIDs.resize(_meshData.nIDs.size() - toerase.size());
	}

	std::vector<std::pair<esint, esint> > buckets;
	if (toinsert.size()) {
		utils::sortAndRemoveDuplicates(toinsert);
		_nregions.resize(_nregsize * (_meshData.nIDs.size() + toinsert.size()));
		_meshData.coordinates.resize(_meshData.nIDs.size() + toinsert.size());
		_meshData.nIDs.resize(_meshData.nIDs.size() + toinsert.size());
		for (size_t i = _meshData.nIDs.size() - toinsert.size(), j = 0; j < toinsert.size(); ++i, ++j) {
			_meshData.nIDs[i] = ids[toinsert[j]];
			_meshData.coordinates[i] = coordinates[toinsert[j]];
			buckets.push_back(std::make_pair(_meshData.nIDs[i], _sfc.getBucket(_meshData.coordinates[i])));
			for (size_t r = 0; r < _nregsize; ++r) {
				_nregions[i * _nregsize + r] = regions[toinsert[j] * _nregsize + r];
			}
		}
	}
	profiler::synccheckpoint("merge");
	{
		// update original buckets
		Communication::allGatherUnknownSize(buckets);
		for (size_t i = 0; i < buckets.size(); ++i) {
			if (_nDistribution[info::mpi::rank] <= buckets[i].first && buckets[i].first < _nDistribution[info::mpi::rank + 1]) {
				auto it = std::lower_bound(_nIDs.begin(), _nIDs.end(), buckets[i].first);
				if (it != _nIDs.end() && *it == buckets[i].first) {
					_nBuckets[it - _nIDs.begin()] = buckets[i].second;
				} else {
					eslog::error("ESPRESO internal error: cannot update duplicated node bucket.\n");
				}
			}
		}
	}

	// we need to sum regions for all duplicate nodes
	// there can be more optimal solution, but this keep the rest of the code unchanged
	searchDuplicateNodes();
	for (auto dup = _meshData._duplicateNodes.begin(); dup != _meshData._duplicateNodes.end(); ++dup) {
		for (size_t i = 0; i < _nregsize; i++) {
			_nregions[_nregsize * dup->toffset + i] |= _nregions[_nregsize * dup->idoffset + i];
		}
	}
	for (auto dup = _meshData._duplicateNodes.begin(); dup != _meshData._duplicateNodes.end(); ++dup) {
		for (size_t i = 0; i < _nregsize; i++) {
			_nregions[_nregsize * dup->idoffset + i] |= _nregions[_nregsize * dup->toffset + i];
		}
	}
	profiler::synccheckpoint("final_search");
	profiler::syncend("merge_duplicated_nodes");
}

void ScatteredInput::linkup()
{
	// 1. Compute neighbors buckets
	// 3. Ask neighbors for coordinates
	// 4. Ask original coordinate holders for the rest nodes (for unknown nodes)
	// 5. Compute nodes neighbors

	// 1. Compute neighbors buckets
//	_sfc.SCFToXYZ();

//	VTKLegacyDebugInfo::spaceFillingCurve(_sfc, _bucketsBorders);

	eslog::startln("LINKUP: CONNECTING CLUSTERS", "LINKUP");
	profiler::syncstart("linkup");

	// 2. Exchange elements having only one node on here

	// 3. Ask neighbors for coordinates

	sortNodes();

	profiler::synccheckpoint("sort");

	// send, found, received
	std::vector<std::vector<esint> > sNodes(_sfcNeighbors.size()), fNodes(_sfcNeighbors.size()), fRegions(_sfcNeighbors.size()), rNodes(_sfcNeighbors.size()), rRegions(_sfcNeighbors.size());
	std::vector<std::vector<Point> > fCoords(_sfcNeighbors.size()), rCoors(_sfcNeighbors.size());

	size_t enodesize = 0;
	size_t estart = info::mesh->dimension == 3 ? 0 : 1;
	for (esint e = 0; e < _etypeDistribution[estart]; e++) {
		enodesize += _meshData.esize[e];
	}

	std::vector<esint> enodes(_meshData.enodes.begin(), _meshData.enodes.begin() + enodesize);
	utils::sortAndRemoveDuplicates(enodes);

	for (size_t id = 0, node = 0; id < _meshData.nIDs.size() || node < enodes.size(); ++id) {
		while (node < enodes.size() && (id == _meshData.nIDs.size() || enodes[node] < _meshData.nIDs[id])) {
			sNodes[0].push_back(enodes[node++]);
		}
		if (node < enodes.size() && enodes[node] == _meshData.nIDs[id]) {
			++node;
		}
	}

	for (size_t t = 1; t < _sfcNeighbors.size(); t++) {
		sNodes[t] = sNodes[0];
	}

	profiler::synccheckpoint("search_unknown");
	eslog::checkpointln("LINKUP: UKNOWN NODES COMPUTED");

	if (!Communication::exchangeUnknownSize(sNodes, rNodes, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: request for coordinates.\n");
	}
	profiler::synccheckpoint("exchange_unknown");
	eslog::checkpointln("LINKUP: UKNOWN NODES EXCHANGED");

	for (size_t t = 0; t < _sfcNeighbors.size(); t++) {
		auto node = _meshData.nIDs.begin();
		for (size_t n = 0; n < rNodes[t].size(); n++) {
			while (node != _meshData.nIDs.end() && *node < rNodes[t][n]) {
				++node;
			}
			if (node != _meshData.nIDs.end() && *node == rNodes[t][n]) {
				fNodes[t].push_back(*node);
				fRegions[t].insert(fRegions[t].end(), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin()), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin() + 1));
				fCoords[t].push_back(_meshData.coordinates[node - _meshData.nIDs.begin()]);
			}
		}
	}

	profiler::synccheckpoint("process_requests");
	eslog::checkpointln("LINKUP: NODES REQUESTS PROCESSED");

	if (!Communication::exchangeUnknownSize(fNodes, rNodes, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: return requested IDs.\n");
	}
	if (!Communication::exchangeUnknownSize(fRegions, rRegions, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: return requested node regions.\n");
	}
	if (!Communication::exchangeUnknownSize(fCoords, rCoors, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: return requested coordinates.\n");
	}
	profiler::synccheckpoint("return_requested");
	eslog::checkpointln("LINKUP: REQUESTED NODES RETURNED");

	// 3.1 Check if all nodes are found

	size_t nodeSize = 0;
	for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
		nodeSize += rNodes[r].size();
	}

	// 4. If there are some unknown nodes we need to ask their original process
	// Here we use _nDistribution and _nIDs that hold old nodes distribution

	// original
	std::vector<int> oTargets, oSources;
	// send, received (that will be asked for coordinates)
	std::vector<std::vector<int> > sTargets, rTargets;
	// unknown
	std::vector<std::vector<esint> > uNodes, uRegions, fIDs, uIDs;
	std::vector<std::vector<Point> > uCoords;

	if (sNodes.size() && nodeSize != sNodes.front().size()) {
		std::vector<esint> found, unknown(sNodes.front().size() - nodeSize);
		for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
			found.insert(found.end(), rNodes[r].begin(), rNodes[r].end());
		}
		utils::sortAndRemoveDuplicates(found);

		std::set_difference(sNodes.front().begin(), sNodes.front().end(), found.begin(), found.end(), unknown.begin());
		sNodes.clear();

		for (size_t i = 0; i < unknown.size(); i++) {
			int trank = std::lower_bound(_nDistribution.begin(), _nDistribution.end(), unknown[i] + 1) - _nDistribution.begin() - 1;
			if (oTargets.size() == 0 || oTargets.back() != trank) {
				oTargets.push_back(trank);
				sNodes.push_back({});
			}
			sNodes.back().push_back(unknown[i]);
		}
		uNodes.resize(sNodes.size());
	} else {
		sNodes.clear();
	}

	profiler::synccheckpoint("process_returned");
	eslog::checkpointln("LINKUP: UNKNOWN NODES PROCESSED");

	if (!Communication::sendVariousTargets(sNodes, uNodes, oTargets, oSources)) {
		eslog::error("ESPRESO internal error: request for unknown nodes.\n");
	}
	profiler::synccheckpoint("ask_for_missing_holders");
	eslog::checkpointln("LINKUP: ASKED FOR MISSING NEIGHBORS");

	sTargets.resize(oSources.size());
	for (size_t t = 0; t < oSources.size(); t++) {
		for (size_t n = 0; n < uNodes[t].size(); n++) {
			auto node = std::lower_bound(_nIDs.begin(), _nIDs.end(), uNodes[t][n]);
			if (node != _nIDs.end() && *node == uNodes[t][n]) {
				sTargets[t].push_back(std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), _nBuckets[node - _nIDs.begin()] + 1) - _bucketsBorders.begin() - 1);
			} else {
				eslog::error("ESPRESO internal error: something wrong happen during link-up phase (request for unknown node).\n");
			}
		}
	}
	profiler::synccheckpoint("search_missing_holders");
	eslog::checkpointln("LINKUP: MISSING NEIGHBORS FOUND");

	if (!Communication::sendVariousTargets(sTargets, rTargets, oSources)) {
		eslog::error("ESPRESO internal error: return requested unknown node targets.\n");
	}
	profiler::synccheckpoint("return_holders");
	eslog::checkpointln("LINKUP: MISSING NEIGHBORS RETURNED");

	uNodes.clear();

	for (size_t i = 1; i < oTargets.size(); i++) {
		sNodes[0].insert(sNodes[0].end(), sNodes[i].begin(), sNodes[i].end());
		rTargets[0].insert(rTargets[0].end(), rTargets[i].begin(), rTargets[i].end());
	}

	if (oTargets.size()) {
		oTargets.clear();
		std::vector<esint> upermutation(sNodes.front().size());
		std::iota(upermutation.begin(), upermutation.end(), 0);
		std::sort(upermutation.begin(), upermutation.end(), [&] (esint i, esint j) {
			if (rTargets[0][i] == rTargets[0][j]) {
				return sNodes[0][i] < sNodes[0][j];
			}
			return rTargets[0][i] < rTargets[0][j];
		});

		for (size_t i = 0; i < upermutation.size(); i++) {
			if (i == 0 || rTargets[0][upermutation[i]] != rTargets[0][upermutation[i - 1]]) {
				oTargets.push_back(rTargets[0][upermutation[i]]);
				uNodes.push_back({});
			}
			uNodes.back().push_back(sNodes[0][upermutation[i]]);
		}
	}

	sNodes.clear();
	oSources.clear();
	sNodes.swap(uNodes);
	profiler::synccheckpoint("missing_prepared");
	eslog::checkpointln("LINKUP: MISSING NEIGHBORS ADDED");

	if (!Communication::sendVariousTargets(sNodes, uNodes, oTargets, oSources)) {
		eslog::error("ESPRESO internal error: request for unknown nodes.\n");
	}
	profiler::synccheckpoint("exchange_missing");
	eslog::checkpointln("LINKUP: ASKED FOR MISSING NODES");

	fCoords.clear();
	fCoords.resize(oSources.size());
	fRegions.clear();
	fRegions.resize(oSources.size());
	for (size_t t = 0; t < oSources.size(); t++) {
		for (size_t n = 0; n < uNodes[t].size(); n++) {
			auto node = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), uNodes[t][n]);
			if (node != _meshData.nIDs.end() && *node == uNodes[t][n]) {
				fCoords[t].push_back(_meshData.coordinates[node - _meshData.nIDs.begin()]);
				fRegions[t].insert(fRegions[t].end(), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin()), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin() + 1));
			} else {
				eslog::error("ESPRESO internal error: something wrong happen during link-up phase.\n");
			}
		}
	}
	profiler::synccheckpoint("add_missing");
	eslog::checkpointln("LINKUP: MISSING NODES COMPUTED");

	if (!Communication::sendVariousTargets(fRegions, uRegions, oSources)) {
		eslog::error("ESPRESO internal error: return requested unknown node regions.\n");
	}
	if (!Communication::sendVariousTargets(fCoords, uCoords, oSources)) {
		eslog::error("ESPRESO internal error: return requested unknown coordinates.\n");
	}
	profiler::synccheckpoint("return_missing");
	eslog::checkpointln("LINKUP: MISSING NODES RETURNED");

	// insert new neighbors to neighbors computed from SFC
	for (size_t i = 0; i < oTargets.size(); i++) {
		auto it = std::lower_bound(_sfcNeighbors.begin(), _sfcNeighbors.end(), oTargets[i]);
		size_t offset = it - _sfcNeighbors.begin();
		_sfcNeighbors.insert(it, oTargets[i]);
		rNodes.insert(rNodes.begin() + offset, sNodes[i]);
		rRegions.insert(rRegions.begin() + offset, uRegions[i]);
		rCoors.insert(rCoors.begin() + offset, uCoords[i]);
		fNodes.insert(fNodes.begin() + offset, std::vector<esint>());
	}

	for (size_t i = 0; i < oSources.size(); i++) {
		auto it = std::lower_bound(_sfcNeighbors.begin(), _sfcNeighbors.end(), oSources[i]);
		size_t offset = it - _sfcNeighbors.begin();
		if (it != _sfcNeighbors.end() && *it == oSources[i]) {
			fNodes[offset].swap(uNodes[i]);
		} else {
			_sfcNeighbors.insert(it, oSources[i]);
			rNodes.insert(rNodes.begin() + offset, std::vector<esint>());
			rRegions.insert(rRegions.begin() + offset, std::vector<esint>());
			rCoors.insert(rCoors.begin() + offset, std::vector<Point>());
			fNodes.insert(fNodes.begin() + offset, uNodes[i]);
		}
	}

	profiler::synccheckpoint("include_missing");
	eslog::checkpointln("LINKUP: MISSING NODES INCLUDED");

	// exchange duplicates nodes and remove them
	if (_meshData.removeDuplicates) {
		std::vector<std::vector<esint> > fDuplicates(_sfcNeighbors.size()), rDuplicates(_sfcNeighbors.size());
		std::vector<esint> dpermutation(_meshData._duplicateNodes.size());
		std::iota(dpermutation.begin(), dpermutation.end(), 0);
		std::sort(dpermutation.begin(), dpermutation.end(), [&] (esint i, esint j) {
			return _meshData._duplicateNodes[i].target < _meshData._duplicateNodes[j].target;
		});

		// project fNodes
		for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
			std::vector<bool> send(_meshData._duplicateNodes.size(), false);
			std::vector<esint> toerase, toinsert;
			auto duplicate = _meshData._duplicateNodes.begin();
			auto tduplicate = dpermutation.begin();
			for (size_t n = 0; n < fNodes[r].size(); n++) {
				auto process = [&] (std::vector<MeshBuilder::Duplicate>::iterator duplicate) {
					send[duplicate - _meshData._duplicateNodes.begin()] = true;
					fDuplicates[r].push_back(duplicate->id);
					fDuplicates[r].push_back(duplicate->target);
					if (std::binary_search(fNodes[r].begin(), fNodes[r].end(), duplicate->id)) {
						toerase.push_back(duplicate->id);
					}
					if (!std::binary_search(fNodes[r].begin(), fNodes[r].end(), duplicate->target)) {
						toinsert.push_back(duplicate->target);
					}
				};

				while (duplicate != _meshData._duplicateNodes.end() && duplicate->id < fNodes[r][n]) {
					++duplicate;
				}
				if (duplicate != _meshData._duplicateNodes.end() && duplicate->id == fNodes[r][n]) {
					if (!send[duplicate - _meshData._duplicateNodes.begin()]) {
						process(duplicate);
						auto pit = std::lower_bound(dpermutation.begin(), dpermutation.end(), duplicate->target, [&] (esint i, esint target) {
							return _meshData._duplicateNodes[i].target < target;
						});
						while (pit != dpermutation.end() && _meshData._duplicateNodes[*pit].target == duplicate->target) {
							if (!send[*pit]) {
								process(_meshData._duplicateNodes.begin() + *pit);
							}
							++pit;
						}
					}
				}
				while (tduplicate != dpermutation.end() &&  _meshData._duplicateNodes[*tduplicate].target < fNodes[r][n]) {
					++tduplicate;
				}
				while (tduplicate != dpermutation.end() &&  _meshData._duplicateNodes[*tduplicate].target == fNodes[r][n]) {
					if (!send[*tduplicate]) {
						process(_meshData._duplicateNodes.begin() + *tduplicate);
					}
					++tduplicate;
				}
			}
			if (toerase.size() || toinsert.size()) {
				utils::sortAndRemoveDuplicates(toerase);
				utils::sortAndRemoveDuplicates(toinsert);
				auto erase = toerase.begin();
				auto insert = toinsert.begin();
				auto node = fNodes[r].begin();
				std::vector<esint> final;
				final.reserve(fNodes[r].size() + toinsert.size() - toerase.size());
				while (final.size() < final.capacity()) {
					if (node == fNodes[r].end()) {
						final.insert(final.end(), insert, toinsert.end());
					} else {
						while (insert != toinsert.end() && *insert < *node) {
							final.push_back(*insert++);
						}
						while (erase != toerase.end() && *erase < *node) { ++erase; }
						if (erase != toerase.end() && *erase == *node) {
							++node; ++erase;
						} else {
							final.push_back(*node++);
						}
					}
				}
				final.swap(fNodes[r]);
			}
		}
		profiler::synccheckpoint("project_duplicates");

		eslog::checkpointln("LINKUP: DUPLICATED NODES PROJECTED");

		if (!Communication::exchangeUnknownSize(fDuplicates, rDuplicates, _sfcNeighbors)) {
			eslog::error("ESPRESO internal error: return duplicate nodes.\n");
		}
		profiler::synccheckpoint("exchange_duplicates");
		eslog::checkpointln("LINKUP: DUPLICATED NODES EXCHANGED");

		for (size_t n = 0; n < _sfcNeighbors.size(); n++) {
			for (size_t i = 0; i < rDuplicates[n].size(); i += 2) {
				_meshData._duplicateNodes.push_back(MeshBuilder::Duplicate{ rDuplicates[n][i], rDuplicates[n][i + 1] });
			}
		}
		std::sort(_meshData._duplicateNodes.begin(), _meshData._duplicateNodes.end(), MeshBuilder::Duplicate());
		profiler::synccheckpoint("insert_duplicates");

		{ // remove duplicates from enodes
			std::vector<esint> toerase, toinsert;
			auto duplicate = _meshData._duplicateNodes.begin();
			for (size_t n = 0; n < enodes.size(); n++) {
				while (duplicate != _meshData._duplicateNodes.end() && duplicate->id < enodes[n]) {
					++duplicate;
				}
				if (duplicate != _meshData._duplicateNodes.end() && duplicate->id == enodes[n]) {
					toerase.push_back(duplicate->id);
					if (!std::binary_search(enodes.begin(), enodes.end(), duplicate->target)) {
						toinsert.push_back(duplicate->target);
					}
				}
			}
			if (toerase.size() || toinsert.size()) {
				utils::sortAndRemoveDuplicates(toerase);
				utils::sortAndRemoveDuplicates(toinsert);
				auto erase = toerase.begin();
				auto insert = toinsert.begin();
				auto node = enodes.begin();
				std::vector<esint> final;
				final.reserve(enodes.size() + toinsert.size() - toerase.size());
				while (final.size() < final.capacity()) {
					if (node == enodes.end()) {
						final.insert(final.end(), insert, toinsert.end());
					} else {
						while (insert != toinsert.end() && *insert < *node) {
							final.push_back(*insert++);
						}
						while (erase != toerase.end() && *erase < *node) { ++erase; }
						if (erase != toerase.end() && *erase == *node) {
							++node; ++erase;
						} else {
							final.push_back(*node++);
						}
					}
				}
				final.swap(enodes);
			}
			profiler::synccheckpoint("remove_duplicates_enodes");
		}

		{ // remove duplicates from IDs
			size_t final = 0;
			auto duplicate = _meshData._duplicateNodes.begin();
			for (size_t n = 0; n < _meshData.nIDs.size(); n++) {
				while (duplicate != _meshData._duplicateNodes.end() && duplicate->id < _meshData.nIDs[n]) {
					++duplicate;
				}
				if (duplicate == _meshData._duplicateNodes.end() || duplicate->id != _meshData.nIDs[n]) {
					for (size_t r = 0; r < _nregsize; ++r) {
						_nregions[final * _nregsize + r] = _nregions[n * _nregsize + r];
					}
					_meshData.nIDs[final] = _meshData.nIDs[n];
					_meshData.coordinates[final] = _meshData.coordinates[n];
					++final;
				}
			}
			_nregions.resize(_nregsize * final);
			_meshData.nIDs.resize(final);
			_meshData.coordinates.resize(final);
			profiler::synccheckpoint("remove_duplicates_ids");
		}

		{ // remove duplicates from returned data
			for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
				auto duplicate = _meshData._duplicateNodes.begin();
				for (size_t n = 0; n < rNodes[r].size(); n++) {
					while (duplicate != _meshData._duplicateNodes.end() && duplicate->id < rNodes[r][n]) {
						++duplicate;
					}
					if (duplicate != _meshData._duplicateNodes.end() && duplicate->id == rNodes[r][n]) {
						auto other = std::lower_bound(rNodes[r].begin(), rNodes[r].end(), duplicate->target);
						if (other != rNodes[r].end() && *other == duplicate->target) {
							// only erase
							rNodes[r].erase(rNodes[r].begin() + n);
							rCoors[r].erase(rCoors[r].begin() + n);
							rRegions[r].erase(rRegions[r].begin() + n * _nregsize, rRegions[r].begin() + n * _nregsize + _nregsize);
							--n;
						} else {
							// move before other (id that is larger than target)
							Point coo = rCoors[r][n];
							std::vector<esint> reg(rRegions[r].begin() + n * _nregsize, rRegions[r].begin() + n * _nregsize + _nregsize);
							size_t upper = other - rNodes[r].begin();
							if (n < upper) {
								for (size_t i = n; i + 1 < upper; ++i) {
									rNodes[r][i] = rNodes[r][i + 1];
									rCoors[r][i] = rCoors[r][i + 1];
									for (size_t j = 0; j < _nregsize; ++j) {
										rRegions[r][_nregsize * i + j] = rRegions[r][_nregsize * (i + 1) + j];
									}
								}
								rNodes[r][upper - 1] = duplicate->target;
								rCoors[r][upper - 1] = coo;
								for (size_t j = 0; j < _nregsize; ++j) {
									rRegions[r][_nregsize * (upper - 1) + j] = reg[j];
								}
								--n;
							} else {
								for (size_t i = n; i > upper; --i) {
									rNodes[r][i] = rNodes[r][i - 1];
									rCoors[r][i] =  rCoors[r][i - 1];
									for (size_t j = 0; j < _nregsize; ++j) {
										rRegions[r][_nregsize * i + j] = rRegions[r][_nregsize * (i - 1) + j];
									}
								}
								rNodes[r][upper] = duplicate->target;
								rCoors[r][upper] = coo;
								for (size_t j = 0; j < _nregsize; ++j) {
									rRegions[r][_nregsize * (upper) + j] = reg[j];
								}
							}
						}
					}
				}
			}
			profiler::synccheckpoint("remove_duplicates_returned");
		}
		coupleDuplicateNodes();

		eslog::checkpointln("LINKUP: DUPLICATED NODES COUPLED");
	}

	_sfcNeighbors.push_back(info::mpi::rank);
	std::sort(_sfcNeighbors.begin(), _sfcNeighbors.end());

	// 5. Compute nodes neighbors
	std::vector<std::vector<esint> > sRanks(_sfcNeighbors.size()), rRanks(_sfcNeighbors.size());

	size_t rankindex;
	std::vector<std::vector<esint> > nodeRequests(_sfcNeighbors.size());
	for (size_t r = 0, i = 0; r < _sfcNeighbors.size(); r++) {
		if (_sfcNeighbors[r] == info::mpi::rank) {
			rankindex = r;
			nodeRequests[r].swap(enodes);
		} else {
			nodeRequests[r].swap(fNodes[i++]);
		}
	}
	std::vector<esint> ranks, ranksOffset;
	std::vector<std::vector<esint>::const_iterator> rPointer(nodeRequests.size());
	for (size_t r = 0; r < nodeRequests.size(); r++) {
		rPointer[r] = nodeRequests[r].begin();
	}
	for (size_t n = 0; n < _meshData.nIDs.size(); ++n) {
		ranks.clear();
		ranksOffset.clear();
		for (size_t r = 0; r < nodeRequests.size(); r++) {
			while (rPointer[r] != nodeRequests[r].end() && *rPointer[r] < _meshData.nIDs[n]) {
				++rPointer[r];
			}
			if (rPointer[r] != nodeRequests[r].end() && *rPointer[r] == _meshData.nIDs[n]) {
				ranksOffset.push_back(r);
				ranks.push_back(_sfcNeighbors[r]);
				++rPointer[r];
			}
		}
		for (size_t r = 0; r < ranks.size(); r++) {
			sRanks[ranksOffset[r]].push_back(ranksOffset.size());
			sRanks[ranksOffset[r]].insert(sRanks[ranksOffset[r]].end(), ranks.begin(), ranks.end());
		}
	}

	nodeRequests[rankindex].swap(enodes);

	// remove nodes without elements
	size_t unique = 0;
	for (size_t id = 0, node = 0; id < _meshData.nIDs.size(); ++id) {
		while (node < enodes.size() && enodes[node] < _meshData.nIDs[id]) {
			++node;
		}
		if (node == enodes.size()) {
			break;
		}
		if (_meshData.nIDs[id] == enodes[node]) {
			_meshData.nIDs[unique] = _meshData.nIDs[id];
			_meshData.coordinates[unique] = _meshData.coordinates[id];
			for (size_t i = 0; i < _nregsize; i++) {
				_nregions[_nregsize * unique + i] = _nregions[_nregsize * id + i];
			}
			++unique;
		}
	}

	_meshData.nIDs.resize(unique);
	_meshData.coordinates.resize(unique);
	_nregions.resize(_nregsize * unique);

	profiler::synccheckpoint("compute_rankmap");
	eslog::checkpointln("LINKUP: NODES RANK MAP COMPUTED");

	if (!Communication::exchangeUnknownSize(sRanks, rRanks, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: exchange ranks data.\n");
	}
	profiler::synccheckpoint("exchange_rankmap");
	eslog::checkpointln("LINKUP: NODES RANK MAP EXCHANGED");

	for (size_t t = 0, i = 0; t < _sfcNeighbors.size(); t++) {
		if (_sfcNeighbors[t] != info::mpi::rank) {
			_meshData.nIDs.insert(_meshData.nIDs.end(), rNodes[i].begin(), rNodes[i].end());
			_nregions.insert(_nregions.end(), rRegions[i].begin(), rRegions[i].end());
			_meshData.coordinates.insert(_meshData.coordinates.end(), rCoors[i].begin(), rCoors[i].end());
			++i;
		}
	}

	size_t r = 0;
	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));

		nregion->second.clear();
		for (size_t i = 0; i < _meshData.nIDs.size(); ++i) {
			if (_nregions[_nregsize * i + byte] & bit) {
				nregion->second.push_back(i);
			}
		}
	}

	_meshData._nrankdist.push_back(0);
	for (size_t n = 0; n < rRanks[rankindex].size(); n += rRanks[rankindex][n] + 1) {
		_meshData._nranks.insert(_meshData._nranks.end(), rRanks[rankindex].begin() + n + 1, rRanks[rankindex].begin() + n + 1 + rRanks[rankindex][n]);
		_meshData._nrankdist.push_back(_meshData._nranks.size());
	}
	for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
		if (_sfcNeighbors[r] != info::mpi::rank) {
			for (size_t n = 0; n < rRanks[r].size(); n += rRanks[r][n] + 1) {
				_meshData._nranks.insert(_meshData._nranks.end(), rRanks[r].begin() + n + 1, rRanks[r].begin() + n + 1 + rRanks[r][n]);
				_meshData._nrankdist.push_back(_meshData._nranks.size());
			}
		}
	}

	profiler::synccheckpoint("merge_rankmaps");
	eslog::endln("LINKUP: LINKED UP");
	profiler::syncend("linkup");
}

void ScatteredInput::exchangeBoundary()
{
	eslog::startln("BOUNDARY: TARGET BOUNDARY ELEMENTS", "BOUNDARY");
	profiler::syncstart("exchange_boundary");

	int threads = info::env::OMP_NUM_THREADS;

	size_t estart = info::mesh->dimension == 3 ? 0 : 1;

	_meshData._edist.clear();
	std::vector<esint> edist = { 0 };
	edist.reserve(_meshData.eIDs.size() - _etypeDistribution[estart] + 1);
	for (esint e = 0; e < _etypeDistribution[estart]; e++) {
		edist.back() += _meshData.esize[e];
	}
	for (esint e = _etypeDistribution[estart]; e < _etypeDistribution.back(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<esint> distribution = tarray<esint>::distribute(threads, _etypeDistribution.back() - _etypeDistribution[estart]);
	std::vector<esint> emembership(distribution.back(), -1);
	std::vector<std::vector<std::pair<esint, esint> > > etargets(threads);
	std::vector<std::vector<esint> > unodes(threads);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		std::vector<esint> nlinks;
		for (esint e = distribution[t]; e < distribution[t + 1]; ++e) {
			nlinks.clear();
			size_t usize = unodes[t].size();
			for (auto n = edist[e]; n < edist[e + 1]; ++n) {
				auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), _meshData.enodes[n]);
				if (nit != _meshData.nIDs.end() && *nit == _meshData.enodes[n]) {
					if (usize == unodes[t].size()) { // nlinks are unused if any node is unknown
						auto links = info::mesh->nodes->elements->cbegin() + (nit - _meshData.nIDs.begin());
						nlinks.insert(nlinks.end(), links->begin(), links->end());
					}
				} else {
					unodes[t].push_back(_meshData.enodes[n]);
				}
			}

			if (usize == unodes[t].size()) { // all nodes are known
				std::sort(nlinks.begin(), nlinks.end());
				int counter = 1;
				for (size_t i = 1; i < nlinks.size(); ++i) {
					if (nlinks[i - 1] == nlinks[i]) {
						++counter;
						if (counter == edist[e + 1] - edist[e]) {
							if (_eDistribution[info::mpi::rank] <= nlinks[i] && nlinks[i] < _eDistribution[info::mpi::rank + 1]) {
								emembership[e] = nlinks[i];
							} else {
								etargets[t].push_back(std::make_pair(std::lower_bound(_eDistribution.begin(), _eDistribution.end(), nlinks[i] + 1) - _eDistribution.begin() - 1, e));
							}
							break;
						}
					} else {
						counter = 1;
					}
				}
			}
		}
	}

	for (int t = 1; t < threads; t++) {
		etargets[0].insert(etargets[0].end(), etargets[t].begin(), etargets[t].end());
		unodes[0].insert(unodes[0].end(), unodes[t].begin(), unodes[t].end());
	}
	utils::sortAndRemoveDuplicates(unodes[0]);
	unodes.resize(_sfcNeighbors.size());
	for (size_t n = 1; n < _sfcNeighbors.size(); n++) {
		unodes[n] = unodes[0];
	}

	/// FIND TARGETS FOR FACES WITH UNKNOWN NODES

	std::vector<std::vector<int> > fLinks(_sfcNeighbors.size()), rLinks(_sfcNeighbors.size());
	std::vector<std::vector<esint> > rnodes(_sfcNeighbors.size());

	profiler::synccheckpoint("compute_unknown_nodes");
	eslog::checkpointln("BOUNDARY: UNKNOWN NODES COMPUTED");

	if (!Communication::exchangeUnknownSize(unodes, rnodes, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: request for unknown boundary nodes.\n");
	}
	profiler::synccheckpoint("exchange_unknown_nodes");

	eslog::checkpointln("BOUNDARY: UNKNOWN NODES EXCHANGED");

	for (size_t t = 0; t < _sfcNeighbors.size(); t++) {
		auto node = _meshData.nIDs.begin();
		for (size_t n = 0; n < rnodes[t].size(); n++) {
			while (node != _meshData.nIDs.end() && *node < rnodes[t][n]) {
				++node;
			}
			if (node != _meshData.nIDs.end() && *node == rnodes[t][n]) {
				auto links = info::mesh->nodes->elements->cbegin() + (node - _meshData.nIDs.begin());
				fLinks[t].push_back(links->size());
				fLinks[t].insert(fLinks[t].end(), links->begin(), links->end());
				fLinks[t].push_back(*node);
			} else {
				if (_meshData.removeDuplicates) {
					auto it = std::lower_bound(_meshData._duplicateNodes.begin(), _meshData._duplicateNodes.end(), rnodes[t][n], [] (MeshBuilder::Duplicate &dup, esint n) {
						return dup.id < n;
					});
					if (it != _meshData._duplicateNodes.end() && it->id == rnodes[t][n]) {
						auto nn = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), it->target);
						if (nn != _meshData.nIDs.end() && *nn == it->target) {
							auto links = info::mesh->nodes->elements->cbegin() + (nn - _meshData.nIDs.begin());
							fLinks[t].push_back(links->size());
							fLinks[t].insert(fLinks[t].end(), links->begin(), links->end());
							fLinks[t].push_back(it->target);
						} else {
							fLinks[t].push_back(-1);
							fLinks[t].push_back(it->target);
						}
					} else {
						fLinks[t].push_back(0);
					}
				} else {
					fLinks[t].push_back(0);
				}
			}
		}
	}

	profiler::synccheckpoint("process_unknown_nodes");
	eslog::checkpointln("BOUNDARY: UNKNOWN NODES PROCESSED");

	if (!Communication::exchangeUnknownSize(fLinks, rLinks, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: return ranks of unknown boundary nodes.\n");
	}
	profiler::synccheckpoint("exchange_links");
	eslog::checkpointln("BOUNDARY: UNKNOWN NODES RETURNED");

	std::vector<MeshBuilder::Duplicate> bduplicate;
	std::vector<esint> found(2 * unodes[0].size(), -1); // [ -1: not found, -2: reindex to my node, >= 0: neighbor node]
	for (size_t n = 0; n < _sfcNeighbors.size(); n++) {
		for (size_t i = 0, noffset = 0; i < found.size(); i += 2) {
			if (rLinks[n][noffset] == -1) {
				if (found[i + 1] == -1) {
					found[i + 1] = rLinks[n][noffset + 1];
					auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), found[i + 1]);
					if (nit != _meshData.nIDs.end() && *nit == found[i + 1]) {
						found[i] = -2;
					}
					bduplicate.push_back(MeshBuilder::Duplicate(unodes[0][i / 2], rLinks[n][noffset + 1]));
				}
				noffset += 2;
				continue;
			}
			if (rLinks[n][noffset] > 0) {
				if (found[i] == -1) {
					found[i] = n;
					found[i + 1] = noffset;
					if (_meshData.removeDuplicates) {
						if (rLinks[n][noffset + rLinks[n][noffset] + 1] != unodes[0][i / 2]) {
							bduplicate.push_back(MeshBuilder::Duplicate(unodes[0][i / 2], rLinks[n][noffset + rLinks[n][noffset] + 1]));
						}
					}
				}
				noffset += rLinks[n][noffset] + 2;
				continue;
			}
			++noffset;
		}
	}

	std::vector<esint> uunodes, rrLinks;
	for (size_t i = 0; i < found.size(); i += 2) {
		if (found[i] == -1) {
			if (found[i + 1] != -1) {
				uunodes.push_back(found[i + 1]);
			} else {
				uunodes.push_back(unodes[0][i / 2]);
			}
		}
	}

	profiler::synccheckpoint("compute_missing_nodes");
	eslog::checkpointln("BOUNDARY: MISSING NODES COMPUTED");

	if (!Communication::allGatherUnknownSize(uunodes)) {
		eslog::error("ESPRESO internal error: allgather unknown nodes.\n");
	}
	profiler::synccheckpoint("allgather_missing_nodes");

	utils::sortAndRemoveDuplicates(uunodes);
	eslog::checkpointln("BOUNDARY: MISSING NODES EXCHAGNED");

	for (size_t i = 0; i < uunodes.size(); i++) {
		auto node = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), uunodes[i]);
		if (node != _meshData.nIDs.end() && *node == uunodes[i]) { // i have the node
			auto ranks = info::mesh->nodes->ranks->begin() + (node - _meshData.nIDs.begin());
			if (ranks->front() == info::mpi::rank) { // i am the first rank that hold the node
				auto links = info::mesh->nodes->elements->cbegin() + (node - _meshData.nIDs.begin());
				rrLinks.push_back(uunodes[i]);
				rrLinks.push_back(uunodes[i]);
				rrLinks.push_back(links->size());
				rrLinks.insert(rrLinks.end(), links->begin(), links->end());
			}
		} else {
			if (_meshData.removeDuplicates) {
				auto it = std::lower_bound(_meshData._duplicateNodes.begin(), _meshData._duplicateNodes.end(), uunodes[i], [] (MeshBuilder::Duplicate &dup, esint n) {
					return dup.id < n;
				});
				if (it != _meshData._duplicateNodes.end() && it->id == uunodes[i]) {
					auto target = std::lower_bound(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().end(), it->target);
					if (target != info::mesh->nodes->IDs->datatarray().end() && *target == it->target) {
						auto ranks = info::mesh->nodes->ranks->begin() + (target - info::mesh->nodes->IDs->datatarray().begin());
						if (ranks->front() == info::mpi::rank) { // i am the first rank that hold the node
							auto links = info::mesh->nodes->elements->cbegin() + (target - info::mesh->nodes->IDs->datatarray().begin());
							rrLinks.push_back(it->id);
							rrLinks.push_back(it->target);
							rrLinks.push_back(links->size());
							rrLinks.insert(rrLinks.end(), links->begin(), links->end());
						}
					}
				}
			}
		}
	}

	profiler::synccheckpoint("found_missing_nodes");
	eslog::checkpointln("BOUNDARY: MISSING NODES FOUND");

	if (!Communication::allGatherUnknownSize(rrLinks)) {
		eslog::error("ESPRESO internal error: allgather unknown nodes links.\n");
	}
	profiler::synccheckpoint("exchange_found_nodes");
	eslog::checkpointln("BOUNDARY: MISSING NODES RETURNED");

	std::vector<esint> permutation(uunodes.size());
	for (size_t i = 0, noffset = 0; i < uunodes.size(); ++i, noffset += rrLinks[noffset + 2] + 3) {
		permutation[i] = noffset;
		if (_meshData.removeDuplicates && rrLinks[noffset] != rrLinks[noffset + 1]) {
			bduplicate.push_back(MeshBuilder::Duplicate(rrLinks[noffset], rrLinks[noffset + 1]));
		}
	}
	if (_meshData.removeDuplicates) {
		std::sort(bduplicate.begin(), bduplicate.end(), MeshBuilder::Duplicate());
	}
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return rrLinks[i] < rrLinks[j]; });

	for (size_t i = 0; i < found.size(); i += 2) {
		if (found[i] == -1) {
			if (found[i + 1] != -1) {
				auto it = std::lower_bound(permutation.begin(), permutation.end(), found[i + 1], [&] (esint i, esint j) {
					return rrLinks[i] < j;
				});
				found[i + 1] = *it + 2;
			} else {
				auto it = std::lower_bound(permutation.begin(), permutation.end(), unodes[0][i / 2], [&] (esint i, esint j) {
					return rrLinks[i] < j;
				});
				found[i + 1] = *it + 2;
			}
		}
	}

	std::vector<esint> linkDist = { 0 }, linkData;
	for (size_t i = 0; i < unodes[0].size(); i++) {
		if (found[2 * i] >= 0) {
			esint lindex = found[2 * i];
			esint loffset = found[2 * i + 1];
			esint lsize = rLinks[lindex][loffset];
			linkDist.push_back(linkDist.back() + lsize);
			linkData.insert(linkData.end(), rLinks[lindex].begin() + loffset + 1, rLinks[lindex].begin() + loffset + 1 + lsize);
		}
		if (found[2 * i] == -1) {
			esint loffset = found[2 * i + 1];
			esint lsize = rrLinks[loffset];
			linkDist.push_back(linkDist.back() + lsize);
			linkData.insert(linkData.end(), rrLinks.begin() + loffset + 1, rrLinks.begin() + loffset + 1 + lsize);
		}
		if (found[2 * i] == -2) {
			auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), found[2 * i + 1]);
			if (nit != _meshData.nIDs.end() && *nit == found[2 * i + 1]) {
				auto links = info::mesh->nodes->elements->cbegin() + (nit - _meshData.nIDs.begin());
				linkDist.push_back(linkDist.back() + links->size());
				linkData.insert(linkData.end(), links->begin(), links->end());
			}
		}
	}
	profiler::synccheckpoint("process_found_nodes");

	{
		std::vector<esint> nlinks;
		for (size_t e = 0; e < emembership.size(); ++e) {
			if (emembership[e] != -1) {
				continue;
			}
			nlinks.clear();
			for (auto n = edist[e]; n < edist[e + 1]; ++n) {
				auto nit = std::lower_bound(unodes[0].begin(), unodes[0].end(), _meshData.enodes[n]);
				if (nit != unodes[0].end() && *nit == _meshData.enodes[n]) {
					nlinks.insert(nlinks.end(), linkData.begin() + linkDist[nit - unodes[0].begin()], linkData.begin() + linkDist[nit - unodes[0].begin() + 1]);
					if (_meshData.removeDuplicates) {
						auto it = std::lower_bound(bduplicate.begin(), bduplicate.end(), _meshData.enodes[n], [] (MeshBuilder::Duplicate &dup, esint n) {
							return dup.id < n;
						});
						if (it != bduplicate.end() && it->id == _meshData.enodes[n]) {
							_meshData.enodes[n] = it->target;
						}
					}
				} else {
					auto it = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), _meshData.enodes[n]);
					if (it != _meshData.nIDs.end() && *it == _meshData.enodes[n]) {
						auto links = info::mesh->nodes->elements->cbegin() + (it - _meshData.nIDs.begin());
						nlinks.insert(nlinks.end(), links->begin(), links->end());
					}
				}
			}
			std::sort(nlinks.begin(), nlinks.end());

			size_t tsize = etargets[0].size();
			int counter = 1;
			for (size_t i = 1; i < nlinks.size(); ++i) {
				if (nlinks[i - 1] == nlinks[i]) {
					++counter;
					if (counter == edist[e + 1] - edist[e]) {
						esint rank = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), nlinks[i] + 1) - _eDistribution.begin() - 1;
						etargets[0].push_back(std::make_pair(rank, e));
						break;
					}
				} else {
					counter = 1;
				}
			}

			if (tsize == etargets[0].size()) {
				eslog::error("MESIO error: parent element not found.\n");
			}
		}
		profiler::synccheckpoint("finish_emembership");
	}

	utils::sortAndRemoveDuplicates(etargets[0]);

	std::vector<int> sRanks;
	std::vector<std::vector<esint> > sBoundary, rBoundary;

	for (size_t e = 0; e < etargets[0].size(); ++e) {
		esint eindex = etargets[0][e].second + _etypeDistribution[estart];
		if (!sRanks.size() || sRanks.back() != etargets[0][e].first) {
			sRanks.push_back(etargets[0][e].first);
			sBoundary.push_back({});
		}
		sBoundary.back().push_back(_meshData.esize[eindex]);
		sBoundary.back().push_back(_meshData.etype[eindex]);
		sBoundary.back().insert(sBoundary.back().end(), _eregions.begin() + eindex * _eregsize, _eregions.begin() + (eindex + 1) * _eregsize);
		sBoundary.back().insert(sBoundary.back().end(), _meshData.enodes.begin() + edist[etargets[0][e].second], _meshData.enodes.begin() + edist[etargets[0][e].second + 1]);
	}

	profiler::synccheckpoint("found_parent_elements");
	eslog::checkpointln("BOUNDARY: PARENT ELEMENTS FOUND");

	if (!Communication::sendVariousTargets(sBoundary, rBoundary, sRanks)) {
		eslog::error("ESPRESO internal error: exchange boundary elements.\n");
	}
	profiler::synccheckpoint("exchange_boundary");
	eslog::checkpointln("BOUNDARY: BOUNDARY EXCHANGED");

	for (size_t r = 1; r < rBoundary.size(); r++) {
		rBoundary[0].insert(rBoundary[0].end(), rBoundary[r].begin(), rBoundary[r].end());
	}

	size_t newsize = 0;
	for (size_t i = 0; rBoundary.size() && i < rBoundary[0].size(); ++newsize) {
		_meshData.esize.push_back(rBoundary[0][i++]);
		edist.push_back(edist.back() + _meshData.esize.back());
		_meshData.etype.push_back(rBoundary[0][i++]);
		_eregions.insert(_eregions.end(), rBoundary[0].begin() + i, rBoundary[0].begin() + i + _eregsize);
		i += _eregsize;
		_meshData.enodes.insert(_meshData.enodes.end(), rBoundary[0].begin() + i, rBoundary[0].begin() + i + _meshData.esize.back());
		i += _meshData.esize.back();
	}

	emembership.resize(emembership.size() + newsize, -1);

	{
		std::vector<esint> nlinks;
		for (size_t e = distribution.back(); e < distribution.back() + newsize; ++e) {
			nlinks.clear();
			for (auto n = edist[e]; n < edist[e + 1]; ++n) {
				auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), _meshData.enodes[n]);
				if (nit != _meshData.nIDs.end() && *nit == _meshData.enodes[n]) {
					auto links = info::mesh->nodes->elements->cbegin() + (nit - _meshData.nIDs.begin());
					nlinks.insert(nlinks.end(), links->begin(), links->end());
				}
			}
			std::sort(nlinks.begin(), nlinks.end());
			int counter = 1;
			for (size_t i = 1; i < nlinks.size(); ++i) {
				if (nlinks[i - 1] == nlinks[i]) {
					++counter;
					if (counter == edist[e + 1] - edist[e]) {
						if (_eDistribution[info::mpi::rank] <= nlinks[i] && nlinks[i] < _eDistribution[info::mpi::rank + 1]) {
							emembership[e] = nlinks[i];
						}
						break;
					}
				} else {
					counter = 1;
				}
			}
		}
	}
	profiler::synccheckpoint("process_boundary");

	std::vector<esint> esize, etype, enodes, ereg;

	for (int i = estart; i < 2; i++) {
		size_t bindex = 0;
		for (esint e = _etypeDistribution[estart]; e < _etypeDistribution.back(); ++e, ++bindex) {
			if (static_cast<int>(Mesh::edata[_meshData.etype[e]].type) == 2 - i && emembership[bindex] != -1) {
				for (auto n = edist[bindex]; n < edist[bindex + 1]; ++n) {
					enodes.push_back(_meshData.enodes[n]);
				}
				esize.push_back(_meshData.esize[e]);
				etype.push_back(_meshData.etype[e]);
				ereg.insert(ereg.end(), _eregions.begin() + _eregsize * e, _eregions.begin() + _eregsize * (e + 1));
			}
		}

		for (size_t e = _etypeDistribution.back(); e < _etypeDistribution.back() + newsize; ++e, ++bindex) {
			if (static_cast<int>(Mesh::edata[_meshData.etype[e]].type) == 2 - i && emembership[bindex] != -1) {
				for (auto n = edist[bindex]; n < edist[bindex + 1]; ++n) {
					enodes.push_back(_meshData.enodes[n]);
				}
				esize.push_back(_meshData.esize[e]);
				etype.push_back(_meshData.etype[e]);
				ereg.insert(ereg.end(), _eregions.begin() + _eregsize * e, _eregions.begin() + _eregsize * (e + 1));
			}
		}
	}

	_meshData.esize.resize(_etypeDistribution[estart]);
	_meshData.etype.resize(_etypeDistribution[estart]);
	_meshData.enodes.resize(edist.front());
	_eregions.resize(_etypeDistribution[estart] * _eregsize);

	_meshData.esize.insert(_meshData.esize.end(), esize.begin(), esize.end());
	_meshData.etype.insert(_meshData.etype.end(), etype.begin(), etype.end());
	_meshData.enodes.insert(_meshData.enodes.end(), enodes.begin(), enodes.end());
	_eregions.insert(_eregions.end(), ereg.begin(), ereg.end());

	_etypeDistribution.clear();
	for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
		_etypeDistribution.push_back(std::lower_bound(_meshData.etype.begin(), _meshData.etype.end(), type, [&] (int e, int type) {
			return static_cast<int>(Mesh::edata[e].type) >= type; }) - _meshData.etype.begin()
		);
	}

	_meshData._edist.clear();

	profiler::synccheckpoint("insert_boundary");
	profiler::syncend("exchange_boundary");
	eslog::endln("BOUNDARY: BOUNDARIES INCLUDED");
}







