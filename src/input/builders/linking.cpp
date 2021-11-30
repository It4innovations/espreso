
#include "builder.utils.h"

#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

#include <algorithm>
#include <numeric>

namespace espreso {
namespace builder {

// some nodes are not at SFC neighboring processes -> fix this situation
void exchangeUnknown(ClusteredMesh &clustered, std::vector<esint> &needed, std::vector<std::vector<esint> > &send, std::vector<std::vector<esint> > &received)
{
	struct __unode__ { // unknown node
		esint offset;
		int rank;

		bool operator<(const __unode__ &other) { return rank == other.rank ? offset < other.offset : rank < other.rank; }
	};

	std::vector<__unode__> unknown;
	{ // compute size of found nodes -> if the size is lower than size of requests, we need to ask for unknown nodes
		size_t size = 0;
		for (size_t r = 0; r < received.size(); ++r) {
			size += received[r].size();
		}
		if (size < needed.size()) {
			std::vector<esint> found;
			found.reserve(size);
			unknown.reserve(needed.size() - size);
			for (size_t r = 0; r < received.size(); ++r) {
				found.insert(found.end(), received[r].begin(), received[r].end());
			}
			std::sort(found.begin(), found.end());
			for (size_t n = 0, m = 0; n < needed.size(); ++n) {
				if (m < found.size() && needed[n] == found[m]) {
					++m;
				} else {
					unknown.push_back({needed[n], info::mpi::rank});
				}
			}
		}
	}

	if (!Communication::allGatherUnknownSize(unknown)) {
		eslog::internalFailure("exchange unknown nodes.\n");
	}
	std::sort(unknown.begin(), unknown.end());
	eslog::checkpointln("LINKUP: UNKNOWN NODES EXCHANGED");

	struct __fnode__: __unode__ {
		int holder;
		_Point<esfloat> coordinate;

		__fnode__() = default;
		__fnode__(const __unode__ &node, int holder, const _Point<esfloat> &coordinate): __unode__(node), holder(holder), coordinate(coordinate) {}

		bool operator<(const __unode__ &other) {
			if (rank == other.rank) {
				return offset < other.offset;
			}
			return rank < other.rank;
		}

		bool operator<(const __fnode__ &other) {
			if (rank == other.rank) {
				if (holder == other.holder) {
					return offset < other.offset;
				}
				return holder < other.holder;
			}
			return rank < other.rank;
		}
	};

	std::vector<__fnode__> found;
	for (size_t i = 0; i < unknown.size(); ++i) {
		auto it = std::lower_bound(clustered.noffsets.begin(), clustered.noffsets.end(), unknown[i].offset);
		if (it != clustered.noffsets.end() && *it == unknown[i].offset) {
			size_t n = clustered.noffsets.begin() - it;
			found.push_back(__fnode__{unknown[i], info::mpi::rank, clustered.coordinates[n]});
			if (clustered.neighbors.back() != unknown[i].rank) {
				clustered.neighbors.push_back(unknown[i].rank);
				send.push_back({ 0 });
			}
			++send.back().front();
			send.back().push_back(unknown[i].offset);
		}
	}
	utils::clearVector(unknown);
	if (found.size()) {
		std::vector<int> permutation(clustered.neighbors.size());
		std::iota(clustered.neighbors.begin(), clustered.neighbors.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return clustered.neighbors[i] < clustered.neighbors[j]; });
		std::vector<std::vector<esint> > _send(send.size());
		std::vector<int> _neighbors;
		for (size_t p = 0; p < permutation.size(); ++p) {
			_send[p].swap(send[permutation[p]]);
			_neighbors[p] = clustered.neighbors[p];
		}
		send.swap(_send);
		clustered.neighbors.swap(_neighbors);
	}

	if (!Communication::allGatherUnknownSize(found)) {
		eslog::internalFailure("exchange found nodes.\n");
	}
	std::sort(found.begin(), found.end());

	{ // insert surely found nodes
		auto begin = std::lower_bound(found.begin(), found.end(), __unode__{0, info::mpi::rank});
		auto end = begin;
		while (end != found.end() && end->rank == info::mpi::rank) { // I have some unknown nodes
			if (clustered.neighbors.back() != end->holder) {
				clustered.neighbors.push_back(end->holder);
				received.push_back({ 0 });
			}
			received.back().push_back(end->offset);
			++end;
		}
		received.back().front() = received.back().size() - 1;
		received.back().resize(received.back().size() + utils::reinterpret_size<esint, _Point<esfloat> >(received.back().front()));
		char *pbuffer = reinterpret_cast<char*>(received.back().data() + received.back().front() + 1);
		for (auto it = begin; it != end; ++it, pbuffer += sizeof(_Point<esfloat>)) {
			memcpy(pbuffer, &it->coordinate, sizeof(_Point<esfloat>));
		}

		if (begin != end) {
			std::vector<int> permutation(clustered.neighbors.size());
			std::iota(clustered.neighbors.begin(), clustered.neighbors.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return clustered.neighbors[i] < clustered.neighbors[j]; });
			std::vector<std::vector<esint> > _rBuffer(received.size());
			std::vector<int> _neighbors;
			for (size_t p = 0; p < permutation.size(); ++p) {
				_rBuffer[p].swap(received[permutation[p]]);
				_neighbors[p] = clustered.neighbors[p];
			}
			received.swap(_rBuffer);
			clustered.neighbors.swap(_neighbors);
		}
	}
	utils::clearVector(found);
	eslog::checkpointln("LINKUP: UNKNOWN NODES INCLUDED");
}

// this method can be optimized by requesting nodes to closest buckets only
void linkup(ClusteredMesh &clustered)
{
	eslog::startln("LINKUP: CONNECTING CLUSTERS", "LINKUP");

	std::vector<esint> needed(clustered.enodes.begin(), clustered.enodes.begin() + clustered.edist[clustered.typeDistribution[clustered.dimension == 3 ? 0 : 1]]);
	std::vector<esint> unknown; // nodes held by other processes
	utils::sortAndRemoveDuplicates(needed);
	for (size_t offset = 0, node = 0; offset < clustered.noffsets.size() || node < needed.size(); ++offset) {
		while (node < needed.size() && (offset == clustered.noffsets.size() || needed[node] < clustered.noffsets[offset])) {
			needed.push_back(needed[node++]);
		}
		if (node < needed.size() && needed[node] == clustered.noffsets[offset]) {
			++node;
		}
	}

	eslog::checkpointln("LINKUP: UKNOWN NODES COMPUTED");

	// probably bottleneck for more than 10k MPI processes -> can be improved by computing sNodes to each process separately
	std::vector<std::vector<esint> > requested(clustered.neighbors.size());
	if (!Communication::exchangeUnknownSize(unknown, requested, clustered.neighbors)) {
		eslog::internalFailure("request for coordinates.\n");
	}
	eslog::checkpointln("LINKUP: UKNOWN NODES EXCHANGED");

	std::vector<std::vector<esint> > send(clustered.neighbors.size()), received(clustered.neighbors.size());
	for (size_t r = 0; r < requested.size(); r++) {
		send[r].push_back(0);
		auto node = clustered.noffsets.begin();
		for (size_t n = 0; n < requested[r].size(); n++) {
			while (node != clustered.noffsets.end() && *node < requested[r][n]) {
				++node;
			}
			if (node != clustered.noffsets.end() && *node == requested[r][n]) {
				send[r].push_back(*node);
			}
		}
		send[r].front() = send[r].size() - 1;
		send[r].resize(send[r].size() + utils::reinterpret_size<esint, _Point<esfloat> >(send[r].front()));
		char *pbuffer = reinterpret_cast<char*>(send[r].data() + send[r].front() + 1);
		for (esint n = 0; n < send[r].front(); ++n, pbuffer += sizeof(_Point<esfloat>)) {
			memcpy(pbuffer, clustered.coordinates.data() + n, sizeof(_Point<esfloat>));
		}
	}
	utils::clearVector(requested);
	eslog::checkpointln("LINKUP: NODES REQUESTS PROCESSED");

	if (!Communication::exchangeUnknownSize(send, received, clustered.neighbors)) {
		eslog::internalFailure("return requested IDs.\n");
	}
	eslog::checkpointln("LINKUP: FOUND NODES RETURNED");


	exchangeUnknown(clustered, needed, send, received);

//	_sfcNeighbors.push_back(info::mpi::rank);
//	std::sort(_sfcNeighbors.begin(), _sfcNeighbors.end());
//
//	// 5. Compute nodes neighbors
//	std::vector<std::vector<esint> > sRanks(_sfcNeighbors.size()), rRanks(_sfcNeighbors.size());
//
//	size_t rankindex;
//	std::vector<std::vector<esint> > nodeRequests(_sfcNeighbors.size());
//	for (size_t r = 0, i = 0; r < _sfcNeighbors.size(); r++) {
//		if (_sfcNeighbors[r] == info::mpi::rank) {
//			rankindex = r;
//			nodeRequests[r].swap(enodes);
//		} else {
//			nodeRequests[r].swap(fNodes[i++]);
//		}
//	}
//	std::vector<esint> ranks, ranksOffset;
//	std::vector<std::vector<esint>::const_iterator> rPointer(nodeRequests.size());
//	for (size_t r = 0; r < nodeRequests.size(); r++) {
//		rPointer[r] = nodeRequests[r].begin();
//	}
//	for (size_t n = 0; n < _meshData.nIDs.size(); ++n) {
//		ranks.clear();
//		ranksOffset.clear();
//		for (size_t r = 0; r < nodeRequests.size(); r++) {
//			while (rPointer[r] != nodeRequests[r].end() && *rPointer[r] < _meshData.nIDs[n]) {
//				++rPointer[r];
//			}
//			if (rPointer[r] != nodeRequests[r].end() && *rPointer[r] == _meshData.nIDs[n]) {
//				ranksOffset.push_back(r);
//				ranks.push_back(_sfcNeighbors[r]);
//				++rPointer[r];
//			}
//		}
//		for (size_t r = 0; r < ranks.size(); r++) {
//			sRanks[ranksOffset[r]].push_back(ranksOffset.size());
//			sRanks[ranksOffset[r]].insert(sRanks[ranksOffset[r]].end(), ranks.begin(), ranks.end());
//		}
//	}
//
//	nodeRequests[rankindex].swap(enodes);
//
//	// remove nodes without elements
//	size_t unique = 0;
//	for (size_t id = 0, node = 0; id < _meshData.nIDs.size(); ++id) {
//		while (node < enodes.size() && enodes[node] < _meshData.nIDs[id]) {
//			++node;
//		}
//		if (node == enodes.size()) {
//			break;
//		}
//		if (_meshData.nIDs[id] == enodes[node]) {
//			_meshData.nIDs[unique] = _meshData.nIDs[id];
//			_meshData.coordinates[unique] = _meshData.coordinates[id];
//			for (size_t i = 0; i < _nregsize; i++) {
//				_nregions[_nregsize * unique + i] = _nregions[_nregsize * id + i];
//			}
//			++unique;
//		}
//	}
//
//	_meshData.nIDs.resize(unique);
//	_meshData.coordinates.resize(unique);
//	_nregions.resize(_nregsize * unique);
//
//	profiler::synccheckpoint("compute_rankmap");
//	eslog::checkpointln("LINKUP: NODES RANK MAP COMPUTED");
//
//	if (!Communication::exchangeUnknownSize(sRanks, rRanks, _sfcNeighbors)) {
//		eslog::internalFailure("exchange ranks data.\n");
//	}
//	profiler::synccheckpoint("exchange_rankmap");
//	eslog::checkpointln("LINKUP: NODES RANK MAP EXCHANGED");
//
//	for (size_t t = 0, i = 0; t < _sfcNeighbors.size(); t++) {
//		if (_sfcNeighbors[t] != info::mpi::rank) {
//			_meshData.nIDs.insert(_meshData.nIDs.end(), rNodes[i].begin(), rNodes[i].end());
//			_nregions.insert(_nregions.end(), rRegions[i].begin(), rRegions[i].end());
//			_meshData.coordinates.insert(_meshData.coordinates.end(), rCoors[i].begin(), rCoors[i].end());
//			++i;
//		}
//	}
//
//	size_t r = 0;
//	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion, ++r) {
//		esint byte = r / (8 * sizeof(esint));
//		esint bit = (esint)1 << (r % (8 * sizeof(esint)));
//
//		nregion->second.clear();
//		for (size_t i = 0; i < _meshData.nIDs.size(); ++i) {
//			if (_nregions[_nregsize * i + byte] & bit) {
//				nregion->second.push_back(i);
//			}
//		}
//	}
//
//	_meshData._nrankdist.push_back(0);
//	for (size_t n = 0; n < rRanks[rankindex].size(); n += rRanks[rankindex][n] + 1) {
//		_meshData._nranks.insert(_meshData._nranks.end(), rRanks[rankindex].begin() + n + 1, rRanks[rankindex].begin() + n + 1 + rRanks[rankindex][n]);
//		_meshData._nrankdist.push_back(_meshData._nranks.size());
//	}
//	for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
//		if (_sfcNeighbors[r] != info::mpi::rank) {
//			for (size_t n = 0; n < rRanks[r].size(); n += rRanks[r][n] + 1) {
//				_meshData._nranks.insert(_meshData._nranks.end(), rRanks[r].begin() + n + 1, rRanks[r].begin() + n + 1 + rRanks[r][n]);
//				_meshData._nrankdist.push_back(_meshData._nranks.size());
//			}
//		}
//	}
//
	eslog::endln("LINKUP: LINKED UP");
}

}
}
