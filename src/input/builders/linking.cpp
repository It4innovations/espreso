
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
static void exchangeUnknown(TemporalMesh<MergedNodes, ClusteredElements> &merged, TemporalMesh<LinkedNodes, ClusteredElements> &linked, std::vector<esint> &adjacent, std::vector<std::vector<esint> > &send, std::vector<std::vector<esint> > &received)
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
		auto it = std::lower_bound(merged.nodes->offsets.begin(), merged.nodes->offsets.end(), unknown[i].offset);
		if (it != merged.nodes->offsets.end() && *it == unknown[i].offset) {
			size_t n = merged.nodes->offsets.begin() - it;
			found.push_back(__fnode__{unknown[i], info::mpi::rank, merged.nodes->coordinates[n]});
			if (linked.nodes->neighbors.back() != unknown[i].rank) {
				linked.nodes->neighbors.push_back(unknown[i].rank);
				send.push_back({ 0 });
			}
			++send.back().front();
			send.back().push_back(unknown[i].offset);
		}
	}
	utils::clearVector(unknown);
	if (found.size()) {
		std::vector<int> permutation(linked.nodes->neighbors.size());
		std::iota(linked.nodes->neighbors.begin(), linked.nodes->neighbors.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return linked.nodes->neighbors[i] < linked.nodes->neighbors[j]; });
		std::vector<std::vector<esint> > _send(send.size());
		std::vector<int> _neighbors;
		for (size_t p = 0; p < permutation.size(); ++p) {
			_send[p].swap(send[permutation[p]]);
			_neighbors[p] = linked.nodes->neighbors[p];
		}
		send.swap(_send);
		linked.nodes->neighbors.swap(_neighbors);
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
//			if (linked.nodes->neighbors.back() != end->holder) {
//				linked.nodes->neighbors.push_back(end->holder);
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
//			std::vector<int> permutation(linked.nodes->neighbors.size());
//			std::iota(linked.nodes->neighbors.begin(), linked.nodes->neighbors.end(), 0);
//			std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return linked.nodes->neighbors[i] < linked.nodes->neighbors[j]; });
//			std::vector<std::vector<esint> > _rBuffer(received.size());
//			std::vector<int> _neighbors;
//			for (size_t p = 0; p < permutation.size(); ++p) {
//				_rBuffer[p].swap(received[permutation[p]]);
//				_neighbors[p] = linked.nodes->neighbors[p];
//			}
//			received.swap(_rBuffer);
//			linked.nodes->neighbors.swap(_neighbors);
//		}
	}
	utils::clearVector(found);
	eslog::checkpointln("LINKUP: UNKNOWN NODES INCLUDED");
}

// this method can be optimized by requesting nodes to closest buckets only
void linkup(TemporalMesh<MergedNodes, ClusteredElements> &merged, TemporalMesh<LinkedNodes, ClusteredElements> &linked)
{
	eslog::startln("LINKUP: CONNECTING CLUSTERS", "LINKUP");

	ivector<esint> required(merged.elements->enodes.begin(), merged.elements->enodes.end());
	std::vector<esint> adjacent; // nodes held by other processes
	utils::sortAndRemoveDuplicates(required);
	for (size_t offset = 0, node = 0; offset < merged.nodes->offsets.size() || node < required.size(); ++offset) {
		while (node < required.size() && (offset == merged.nodes->offsets.size() || required[node] < merged.nodes->offsets[offset])) {
			adjacent.push_back(required[node++]);
		}
		if (node < required.size() && required[node] == merged.nodes->offsets[offset]) {
			++node;
		}
	}

	eslog::checkpointln("LINKUP: ADJACENT NODES COMPUTED");

	// probably bottleneck for more than 10k MPI processes -> can be improved by computing sNodes to each process separately
	std::vector<std::vector<esint> > requested(linked.nodes->neighbors.size());
	if (!Communication::exchangeUnknownSize(adjacent, requested, linked.nodes->neighbors)) {
		eslog::internalFailure("request for coordinates.\n");
	}
	eslog::checkpointln("LINKUP: ADJACENT NODES EXCHANGED");

	std::vector<std::vector<esint> > send(linked.nodes->neighbors.size()), received(linked.nodes->neighbors.size());
	for (size_t r = 0; r < requested.size(); r++) {
		send[r].push_back(0);
		for (size_t n = 0; n < requested[r].size(); n++) {
			auto nit = std::lower_bound(merged.nodes->offsets.begin(), merged.nodes->offsets.end(), requested[r][n]);
			if (nit != merged.nodes->offsets.end() && *nit == requested[r][n]) {
				send[r].push_back(nit - merged.nodes->offsets.begin());
			}
		}
		send[r].front() = send[r].size() - 1;
		send[r].resize(send[r].size() + utils::reinterpret_size<esint, _Point<esfloat> >(send[r].front()));
		char *pbuffer = reinterpret_cast<char*>(send[r].data() + send[r].front() + 1);
		for (esint n = 0; n < send[r].front(); ++n, pbuffer += sizeof(_Point<esfloat>)) {
			memcpy(pbuffer, merged.nodes->coordinates.data() + send[r][n + 1], sizeof(_Point<esfloat>));
			send[r][n + 1] = merged.nodes->offsets[send[r][n + 1]];
		}
	}
	eslog::checkpointln("LINKUP: NODES REQUESTS PROCESSED");

	if (!Communication::exchangeUnknownSize(send, received, linked.nodes->neighbors)) {
		eslog::internalFailure("return requested IDs.\n");
	}
	eslog::checkpointln("LINKUP: ADJACENT NODES RETURNED");

	exchangeUnknown(merged, linked, adjacent, send, received);

	std::vector<esint> rankDistribution(merged.nodes->offsets.size() + 1), duplicationIndex(merged.nodes->offsets.size());
	for (size_t i = 0, j = 0; i < required.size(); ++i) { // required size is proportional to clustered.noffsets
		while (merged.nodes->offsets[j] < required[i]) { ++j; }
		if (merged.nodes->offsets[j] == required[i]) {
			++rankDistribution[j];
		}
	}
	for (size_t r = 0; r < send.size(); ++r) {
		for (esint i = 1; i <= send[r].front(); ++i) { // send[r] is much smaller than clustered.noffsets
			size_t n = std::lower_bound(merged.nodes->offsets.begin(), merged.nodes->offsets.end(), send[r][i]) - merged.nodes->offsets.begin();
			++rankDistribution[n];
		}
	}
	for (size_t i = 0; i < merged.nodes->duplication.size(); ++i) {
		// real number of ranks can be lower due to the same ranks in different duplicated nodes
		rankDistribution[merged.nodes->duplication[i].origin] += rankDistribution[merged.nodes->duplication[i].duplication];
		duplicationIndex[merged.nodes->duplication[i].origin] = i;      // always last occurrence
		duplicationIndex[merged.nodes->duplication[i].duplication] = i; // always pointer to offset
	}
	for (size_t i = 0; i < merged.nodes->duplication.size(); ++i) {
		rankDistribution[merged.nodes->duplication[i].duplication] = rankDistribution[merged.nodes->duplication[i].origin];
	}
	utils::sizesToOffsets(rankDistribution);
	std::vector<int, initless_allocator<int> > rankData(rankDistribution.back());

	{ // put all ranks to the rankData
		auto _rankDistribution = rankDistribution;
		size_t n = 0;
		for ( ; n < linked.nodes->neighbors.size() && linked.nodes->neighbors[n] < info::mpi::rank; ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				size_t ni = std::lower_bound(merged.nodes->offsets.begin(), merged.nodes->offsets.end(), send[n][i]) - merged.nodes->offsets.begin();
				rankData[_rankDistribution[ni]++] = linked.nodes->neighbors[n];
			}
		}
		for (size_t i = 0, j = 0; i < required.size(); ++i) { // required size is proportional to clustered.noffsets
			while (merged.nodes->offsets[j] < required[i]) { ++j; }
			if (merged.nodes->offsets[j] == required[i]) {
				rankData[_rankDistribution[j]++] = info::mpi::rank;
			}
		}
		for ( ; n < linked.nodes->neighbors.size(); ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				size_t ni = std::lower_bound(merged.nodes->offsets.begin(), merged.nodes->offsets.end(), send[n][i]) - merged.nodes->offsets.begin();
				rankData[_rankDistribution[ni]++] = linked.nodes->neighbors[n];
			}
		}

		for (size_t i = 0; i < merged.nodes->duplication.size(); ++i) {
			const DataDuplication &dd = merged.nodes->duplication[i];
			for (esint k = rankDistribution[dd.origin], j = rankDistribution[dd.duplication]; j < _rankDistribution[dd.duplication]; ++j) {
				while (k < _rankDistribution[dd.origin] && rankData[k] < rankData[j]) { ++k; }
				if (k == _rankDistribution[dd.origin] || rankData[k] != rankData[j]) {
					rankData[_rankDistribution[dd.origin]++] = rankData[j];
				}
			}
			std::sort(rankData.begin() + rankDistribution[dd.origin], rankData.begin() + _rankDistribution[dd.origin]);
		}
		for (size_t i = 0; i < merged.nodes->duplication.size(); ++i) {
			const DataDuplication &dd = merged.nodes->duplication[i];
			for (esint j = rankDistribution[dd.origin], k = rankDistribution[dd.duplication]; j < _rankDistribution[dd.origin]; ++j, ++k) {
				rankData[k] = rankData[j];
			}
			_rankDistribution[dd.duplication] = rankDistribution[dd.duplication] + _rankDistribution[dd.origin] - rankDistribution[dd.origin];
		}

		size_t sum = 0;
		for (size_t i = 0, j = 0; i < merged.nodes->offsets.size(); ++i, sum = j) {
			for (esint k = rankDistribution[i]; k < _rankDistribution[i]; ++k) {
				rankData[j++] = rankData[k];
			}
			rankDistribution[i] = sum;
		}
		rankDistribution.back() = sum;
		rankData.resize(rankDistribution.back());
	}

	std::vector<std::vector<esint> > ranks(linked.nodes->neighbors.size());
	{ // send rank data and duplication info
		std::vector<std::vector<esint> > sBuffer(linked.nodes->neighbors.size());

		for (size_t n = 0; n < send.size(); ++n) {
			size_t size = 0;
			for (esint i = 1; i <= send[n].front(); ++i) {
				size += 2; // ranks, duplications
				esint ni = std::lower_bound(merged.nodes->offsets.begin(), merged.nodes->offsets.end(), send[n][i]) - merged.nodes->offsets.begin();
				size += rankDistribution[ni + 1] - rankDistribution[ni];
				if (duplicationIndex[ni] != -1) {
					if (ni == merged.nodes->duplication[duplicationIndex[ni]].duplication) {
						size += 2;
					} else {
						esint di = duplicationIndex[ni];
						while (di >= 0 && ni == merged.nodes->duplication[di].origin) {
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
				esint ni = std::lower_bound(merged.nodes->offsets.begin(), merged.nodes->offsets.end(), send[n][i]) - merged.nodes->offsets.begin();
				size_t prevsize = sBuffer[n].size();
				sBuffer[n].push_back(rankDistribution[ni + 1] - rankDistribution[ni]);
				sBuffer[n].push_back(0);
				for (esint r = rankDistribution[ni]; r < rankDistribution[ni + 1]; ++r) {
					sBuffer[n].push_back(rankData[r]);
				}
				if (duplicationIndex[ni] != -1) {
					if (ni == merged.nodes->duplication[duplicationIndex[ni]].duplication) {
						++sBuffer[n][prevsize + 1];
						sBuffer[n].push_back(merged.nodes->offsets[merged.nodes->duplication[duplicationIndex[ni]].origin]);
						sBuffer[n].push_back(merged.nodes->offsets[merged.nodes->duplication[duplicationIndex[ni]].duplication]);
					} else {
						esint di = duplicationIndex[ni];
						while (di >= 0 && ni == merged.nodes->duplication[di].origin) {
							++sBuffer[n][prevsize + 1];
							sBuffer[n].push_back(merged.nodes->offsets[merged.nodes->duplication[di].origin]);
							sBuffer[n].push_back(merged.nodes->offsets[merged.nodes->duplication[di].duplication]);
							--di;
						}
					}
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sBuffer, ranks, linked.nodes->neighbors)) {
			eslog::internalFailure("cannot exchange clustered nodes ranks.\n");
		}
	}

	{ // build nodes
		linked.nodes->offsets.swap(required);
		linked.nodes->coordinates.reserve(linked.nodes->offsets.size());
		linked.nodes->rankDistribution.reserve(linked.nodes->offsets.size() + 1);
		linked.nodes->rankDistribution.push_back(0);
		for (auto dup = merged.nodes->duplication.begin(); dup != merged.nodes->duplication.end(); ++dup) {
			linked.nodes->duplication.push_back({ merged.nodes->offsets[dup->origin], merged.nodes->offsets[dup->duplication] }); // it can include unnecessary nodes
		}

		std::vector<esint> recvOffset(linked.nodes->neighbors.size()), rankOffset(linked.nodes->neighbors.size());
		for (size_t i = 0, j = 0; i < linked.nodes->offsets.size(); ++i) {
			while (merged.nodes->offsets[j] < linked.nodes->offsets[i]) { ++j; }
			if (merged.nodes->offsets[j] == linked.nodes->offsets[i]) {
				linked.nodes->coordinates.push_back(merged.nodes->coordinates[j]);
				for (esint n = rankDistribution[j]; n < rankDistribution[j + 1]; ++n) {
					linked.nodes->rankData.push_back(rankData[n]);
				}
			} else {
				for (size_t r = 0; r < recvOffset.size(); ++r) {
					if (recvOffset[r] < received[r].front() && linked.nodes->offsets[i] == received[r][recvOffset[r] + 1]) {
						_Point<esfloat> *coo = reinterpret_cast<_Point<esfloat>*>(received[r].data() + received[r].front() + 1);
						linked.nodes->coordinates.push_back(coo[recvOffset[r]]);
						++recvOffset[r];
						esint nranks = ranks[r][rankOffset[r]++];
						esint nduplications = ranks[r][rankOffset[r]++];
						for (esint n = 0; n < nranks; ++n) {
							linked.nodes->rankData.push_back(ranks[r][rankOffset[r]++]);
						}
						for (esint n = 0; n < nduplications; ++n, rankOffset[r] += 2) {
							linked.nodes->duplication.push_back({ ranks[r][rankOffset[r]], ranks[r][rankOffset[r] + 1] });
						}
					}
				}
			}
			linked.nodes->rankDistribution.push_back(linked.nodes->rankData.size());
		}
		utils::sortAndRemoveDuplicates(linked.nodes->duplication);

		linked.nodes->g2l.reserve(linked.nodes->offsets.size() + linked.nodes->duplication.size());
		for (size_t n = 0; n < linked.nodes->offsets.size(); ++n) {
			linked.nodes->g2l[linked.nodes->offsets[n]] = n;
		}
		for (auto begin = linked.nodes->duplication.begin(), end = begin; end != linked.nodes->duplication.end(); begin = end) {
			esint known = begin->duplication;
			while (end != linked.nodes->duplication.end() && begin->origin == end->origin) {
				if (known == begin->duplication) {
					if (begin->origin == end->origin && linked.nodes->g2l.find(end->duplication) != linked.nodes->g2l.end()) {
						known = end->duplication;
					}
					if (begin->origin == end->origin && linked.nodes->g2l.find(end->origin) != linked.nodes->g2l.end()) {
						known = end->origin;
					}
				}
				++end;
			}
			for (auto it = begin; it != end; ++it) {
				linked.nodes->g2l[it->origin] = linked.nodes->g2l[known];
				linked.nodes->g2l[it->duplication] = linked.nodes->g2l[known];
			}
		}
	}

	{ // build elements
		linked.elements->offsets.swap(merged.elements->offsets);
		linked.elements->etype.swap(merged.elements->etype);
		linked.elements->enodes.swap(merged.elements->enodes);
		linked.elements->edist.swap(merged.elements->edist);
	}

	merged.clear();
	eslog::endln("LINKUP: LINKED UP");
}

}
}
