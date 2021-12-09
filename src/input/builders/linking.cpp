
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
static void exchangeUnknown(ClusteredMesh &clustered, std::vector<esint> &adjacent, std::vector<std::vector<esint> > &send, std::vector<std::vector<esint> > &received)
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
			size += received[r].front();
		}
		if (size < adjacent.size()) {
			std::vector<esint> found;
			found.reserve(size);
			unknown.reserve(adjacent.size() - size);
			for (size_t r = 0; r < received.size(); ++r) {
				found.insert(found.end(), received[r].begin() + 1, received[r].begin() + 1 + received[r].front());
			}
			std::sort(found.begin(), found.end());
			for (size_t n = 0, m = 0; n < adjacent.size(); ++n) {
				if (m < found.size() && adjacent[n] == found[m]) {
					++m;
				} else {
					unknown.push_back({adjacent[n], info::mpi::rank});
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

	if (found.size()) { // insert surely found nodes
		eslog::internalFailure("insert surely found nodes\n");
//		auto begin = std::lower_bound(found.begin(), found.end(), __unode__{0, info::mpi::rank});
//		auto end = begin;
//		while (end != found.end() && end->rank == info::mpi::rank) { // I have some unknown nodes
//			if (clustered.neighbors.back() != end->holder) {
//				clustered.neighbors.push_back(end->holder);
//				received.push_back({ 0 });
//			}
//			received.back().push_back(end->offset);
//			++end;
//		}
//		received.back().front() = received.back().size() - 1;
//		received.back().resize(received.back().size() + utils::reinterpret_size<esint, _Point<esfloat> >(received.back().front()));
//		char *pbuffer = reinterpret_cast<char*>(received.back().data() + received.back().front() + 1);
//		for (auto it = begin; it != end; ++it, pbuffer += sizeof(_Point<esfloat>)) {
//			memcpy(pbuffer, &it->coordinate, sizeof(_Point<esfloat>));
//		}
//
//		if (begin != end) {
//			std::vector<int> permutation(clustered.neighbors.size());
//			std::iota(clustered.neighbors.begin(), clustered.neighbors.end(), 0);
//			std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return clustered.neighbors[i] < clustered.neighbors[j]; });
//			std::vector<std::vector<esint> > _rBuffer(received.size());
//			std::vector<int> _neighbors;
//			for (size_t p = 0; p < permutation.size(); ++p) {
//				_rBuffer[p].swap(received[permutation[p]]);
//				_neighbors[p] = clustered.neighbors[p];
//			}
//			received.swap(_rBuffer);
//			clustered.neighbors.swap(_neighbors);
//		}
	}
	utils::clearVector(found);
	eslog::checkpointln("LINKUP: UNKNOWN NODES INCLUDED");
}

// this method can be optimized by requesting nodes to closest buckets only
void linkup(ClusteredMesh &clustered, ClusteredMesh &linked)
{
	eslog::startln("LINKUP: CONNECTING CLUSTERS", "LINKUP");

	std::vector<esint> required(clustered.enodes.begin(), clustered.enodes.end());
	std::vector<esint> adjacent; // nodes held by other processes
	utils::sortAndRemoveDuplicates(required);
	for (size_t offset = 0, node = 0; offset < clustered.noffsets.size() || node < required.size(); ++offset) {
		while (node < required.size() && (offset == clustered.noffsets.size() || required[node] < clustered.noffsets[offset])) {
			adjacent.push_back(required[node++]);
		}
		if (node < required.size() && required[node] == clustered.noffsets[offset]) {
			++node;
		}
	}

	eslog::checkpointln("LINKUP: ADJACENT NODES COMPUTED");

	// probably bottleneck for more than 10k MPI processes -> can be improved by computing sNodes to each process separately
	std::vector<std::vector<esint> > requested(clustered.neighbors.size());
	if (!Communication::exchangeUnknownSize(adjacent, requested, clustered.neighbors)) {
		eslog::internalFailure("request for coordinates.\n");
	}
	eslog::checkpointln("LINKUP: ADJACENT NODES EXCHANGED");

	std::vector<std::vector<esint> > send(clustered.neighbors.size()), received(clustered.neighbors.size());
	for (size_t r = 0; r < requested.size(); r++) {
		send[r].push_back(0);
		for (size_t n = 0; n < requested[r].size(); n++) {
			auto nit = std::lower_bound(clustered.noffsets.begin(), clustered.noffsets.end(), requested[r][n]);
			if (nit != clustered.noffsets.end() && *nit == requested[r][n]) {
				send[r].push_back(nit - clustered.noffsets.begin());
			}
		}
		send[r].front() = send[r].size() - 1;
		send[r].resize(send[r].size() + utils::reinterpret_size<esint, _Point<esfloat> >(send[r].front()));
		char *pbuffer = reinterpret_cast<char*>(send[r].data() + send[r].front() + 1);
		for (esint n = 0; n < send[r].front(); ++n, pbuffer += sizeof(_Point<esfloat>)) {
			memcpy(pbuffer, clustered.coordinates.data() + send[r][n + 1], sizeof(_Point<esfloat>));
			send[r][n + 1] = clustered.noffsets[send[r][n + 1]];
		}
	}
	eslog::checkpointln("LINKUP: NODES REQUESTS PROCESSED");

	if (!Communication::exchangeUnknownSize(send, received, clustered.neighbors)) {
		eslog::internalFailure("return requested IDs.\n");
	}
	eslog::checkpointln("LINKUP: ADJACENT NODES RETURNED");

	exchangeUnknown(clustered, adjacent, send, received);

	std::vector<esint> rankDistribution(clustered.noffsets.size() + 1), duplicationIndex(clustered.noffsets.size());
	for (size_t i = 0, j = 0; i < required.size(); ++i) { // required size is proportional to clustered.noffsets
		while (clustered.noffsets[j] < required[i]) { ++j; }
		if (clustered.noffsets[j] == required[i]) {
			++rankDistribution[j];
		}
	}
	for (size_t r = 0; r < send.size(); ++r) {
		for (esint i = 1; i <= send[r].front(); ++i) { // send[r] is much smaller than clustered.noffsets
			size_t n = std::lower_bound(clustered.noffsets.begin(), clustered.noffsets.end(), send[r][i]) - clustered.noffsets.begin();
			++rankDistribution[n];
		}
	}
	for (size_t i = 0; i < clustered.nduplication.size(); ++i) {
		// real number of ranks can be lower due to the same ranks in different duplicated nodes
		rankDistribution[clustered.nduplication[i].origin] += rankDistribution[clustered.nduplication[i].duplication];
		duplicationIndex[clustered.nduplication[i].origin] = i;      // always last occurrence
		duplicationIndex[clustered.nduplication[i].duplication] = i; // always pointer to offset
	}
	for (size_t i = 0; i < clustered.nduplication.size(); ++i) {
		rankDistribution[clustered.nduplication[i].duplication] = rankDistribution[clustered.nduplication[i].origin];
	}
	utils::sizesToOffsets(rankDistribution);
	std::vector<int, initless_allocator<int> > rankData(rankDistribution.back());

	{ // put all ranks to the rankData
		auto _rankDistribution = rankDistribution;
		size_t n = 0;
		for ( ; n < clustered.neighbors.size() && clustered.neighbors[n] < info::mpi::rank; ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				size_t ni = std::lower_bound(clustered.noffsets.begin(), clustered.noffsets.end(), send[n][i]) - clustered.noffsets.begin();
				rankData[_rankDistribution[ni]++] = clustered.neighbors[n];
			}
		}
		for (size_t i = 0, j = 0; i < required.size(); ++i) { // required size is proportional to clustered.noffsets
			while (clustered.noffsets[j] < required[i]) { ++j; }
			if (clustered.noffsets[j] == required[i]) {
				rankData[_rankDistribution[j]++] = info::mpi::rank;
			}
		}
		for ( ; n < clustered.neighbors.size(); ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				size_t ni = std::lower_bound(clustered.noffsets.begin(), clustered.noffsets.end(), send[n][i]) - clustered.noffsets.begin();
				rankData[_rankDistribution[ni]++] = clustered.neighbors[n];
			}
		}

		for (size_t i = 0; i < clustered.nduplication.size(); ++i) {
			const DataDuplication &dd = clustered.nduplication[i];
			for (esint k = rankDistribution[dd.origin], j = rankDistribution[dd.duplication]; j < _rankDistribution[dd.duplication]; ++j) {
				while (k < _rankDistribution[dd.origin] && rankData[k] < rankData[j]) { ++k; }
				if (k == _rankDistribution[dd.origin] || rankData[k] != rankData[j]) {
					rankData[_rankDistribution[dd.origin]++] = rankData[j];
				}
			}
			std::sort(rankData.begin() + rankDistribution[dd.origin], rankData.begin() + _rankDistribution[dd.origin]);
		}
		for (size_t i = 0; i < clustered.nduplication.size(); ++i) {
			const DataDuplication &dd = clustered.nduplication[i];
			for (esint j = rankDistribution[dd.origin], k = rankDistribution[dd.duplication]; j < _rankDistribution[dd.origin]; ++j, ++k) {
				rankData[k] = rankData[j];
			}
			_rankDistribution[dd.duplication] = rankDistribution[dd.duplication] + _rankDistribution[dd.origin] - rankDistribution[dd.origin];
		}

		size_t sum = 0;
		for (size_t i = 0, j = 0; i < clustered.noffsets.size(); ++i, sum = j) {
			for (esint k = rankDistribution[i]; k < _rankDistribution[i]; ++k) {
				rankData[j++] = rankData[k];
			}
			rankDistribution[i] = sum;
		}
		rankDistribution.back() = sum;
		rankData.resize(rankDistribution.back());
	}

	std::vector<std::vector<esint> > ranks(clustered.neighbors.size());
	{ // send rank data and duplication info
		std::vector<std::vector<esint> > sBuffer(clustered.neighbors.size());

		for (size_t n = 0; n < send.size(); ++n) {
			size_t size = 0;
			for (esint i = 1; i <= send[n].front(); ++i) {
				size += 2; // ranks, duplications
				esint ni = std::lower_bound(clustered.noffsets.begin(), clustered.noffsets.end(), send[n][i]) - clustered.noffsets.begin();
				size += rankDistribution[ni + 1] - rankDistribution[ni];
				if (duplicationIndex[ni] != -1) {
					if (ni == clustered.nduplication[duplicationIndex[ni]].duplication) {
						size += 2;
					} else {
						esint di = duplicationIndex[ni];
						while (di >= 0 && ni == clustered.nduplication[di].origin) {
							size += 2;
							--di;
						}
					}
				}
			}
			sBuffer[n].reserve(size);
		}
		for (size_t n = 0; n < send.size(); ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				esint ni = std::lower_bound(clustered.noffsets.begin(), clustered.noffsets.end(), send[n][i]) - clustered.noffsets.begin();
				size_t prevsize = sBuffer[n].size();
				sBuffer[n].push_back(rankDistribution[ni + 1] - rankDistribution[ni]);
				sBuffer[n].push_back(0);
				for (esint r = rankDistribution[ni]; r < rankDistribution[ni + 1]; ++r) {
					sBuffer[n].push_back(rankData[r]);
				}
				if (duplicationIndex[ni] != -1) {
					if (ni == clustered.nduplication[duplicationIndex[ni]].duplication) {
						++sBuffer[n][prevsize + 1];
						sBuffer[n].push_back(clustered.noffsets[clustered.nduplication[duplicationIndex[ni]].origin]);
						sBuffer[n].push_back(clustered.noffsets[clustered.nduplication[duplicationIndex[ni]].duplication]);
					} else {
						esint di = duplicationIndex[ni];
						while (di >= 0 && ni == clustered.nduplication[di].origin) {
							++sBuffer[n][prevsize + 1];
							sBuffer[n].push_back(clustered.noffsets[clustered.nduplication[di].origin]);
							sBuffer[n].push_back(clustered.noffsets[clustered.nduplication[di].duplication]);
							--di;
						}
					}
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sBuffer, ranks, clustered.neighbors)) {
			eslog::internalFailure("cannot exchange clustered nodes ranks.\n");
		}
	}

	linked.dimension = clustered.dimension;
	linked.neighbors = clustered.neighbors;

	{ // build nodes
		linked.noffsets.swap(required);
		linked.coordinates.reserve(linked.noffsets.size());
		linked.rankDistribution.reserve(linked.noffsets.size() + 1);
		linked.rankDistribution.push_back(0);
		for (auto dup = clustered.nduplication.begin(); dup != clustered.nduplication.end(); ++dup) {
			linked.nduplication.push_back({ clustered.noffsets[dup->origin], clustered.noffsets[dup->duplication] }); // it can include unnecessary nodes
		}

		std::vector<esint> recvOffset(clustered.neighbors.size()), rankOffset(clustered.neighbors.size());
		for (size_t i = 0, j = 0; i < linked.noffsets.size(); ++i) {
			while (clustered.noffsets[j] < linked.noffsets[i]) { ++j; }
			if (clustered.noffsets[j] == linked.noffsets[i]) {
				linked.coordinates.push_back(clustered.coordinates[j]);
				for (esint n = rankDistribution[j]; n < rankDistribution[j + 1]; ++n) {
					linked.rankData.push_back(rankData[n]);
				}
			} else {
				for (size_t r = 0; r < recvOffset.size(); ++r) {
					if (recvOffset[r] < received[r].front() && linked.noffsets[i] == received[r][recvOffset[r] + 1]) {
						_Point<esfloat> *coo = reinterpret_cast<_Point<esfloat>*>(received[r].data() + received[r].front() + 1);
						linked.coordinates.push_back(coo[recvOffset[r]]);
						++recvOffset[r];
						esint nranks = ranks[r][rankOffset[r]++];
						esint nduplications = ranks[r][rankOffset[r]++];
						for (esint n = 0; n < nranks; ++n) {
							linked.rankData.push_back(ranks[r][rankOffset[r]++]);
						}
						for (esint n = 0; n < nduplications; ++n, rankOffset[r] += 2) {
							linked.nduplication.push_back({ ranks[r][rankOffset[r]], ranks[r][rankOffset[r] + 1] });
						}
					}
				}
			}
			linked.rankDistribution.push_back(linked.rankData.size());
		}
		std::sort(linked.nduplication.begin(), linked.nduplication.end());
		utils::clearVector(clustered.noffsets);

		linked.g2l.reserve(linked.noffsets.size() + linked.nduplication.size());
		for (size_t n = 0; n < linked.noffsets.size(); ++n) {
			linked.g2l[linked.noffsets[n]] = n;
		}
		for (auto begin = linked.nduplication.begin(), end = begin; end != linked.nduplication.end(); begin = end) {
			esint known = begin->duplication;
			while (end != linked.nduplication.end() && begin->origin == end->origin) {
				if (known == begin->duplication) {
					if (begin->origin == end->origin && linked.g2l.find(end->duplication) != linked.g2l.end()) {
						known = end->duplication;
					}
					if (begin->origin == end->origin && linked.g2l.find(end->origin) != linked.g2l.end()) {
						known = end->origin;
					}
				}
				++end;
			}
			for (auto it = begin; it != end; ++it) {
				linked.g2l[it->origin] = linked.g2l[known];
				linked.g2l[it->duplication] = linked.g2l[known];
			}
		}
	}

	{ // build elements
		linked.eoffsets.swap(clustered.eoffsets);
		linked.etype.swap(clustered.etype);
		linked.enodes.swap(clustered.enodes);
		linked.edist.swap(clustered.edist);
	}
	eslog::endln("LINKUP: LINKED UP");
}

}
}
