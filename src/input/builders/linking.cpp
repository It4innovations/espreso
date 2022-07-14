
#include "builder.utils.h"

#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <queue>

namespace espreso {
namespace builder {

static void exchangeUnknown(MergedNodes &merged, LinkedNodes &linked, std::vector<esint> &adjacent, std::vector<std::vector<esint> > &send, std::vector<std::vector<esint> > &recv)
{
	auto permute = [&] () {
		ivector<int> permutation(linked.neighbors.size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return linked.neighbors[i] < linked.neighbors[j]; });
		std::vector<std::vector<esint> > _send(permutation.size()), _recv(permutation.size());
		std::vector<int> _neighbors(permutation.size());
		for (size_t p = 0; p < permutation.size(); ++p) {
			_send[p].swap(send[permutation[p]]);
			_recv[p].swap(recv[permutation[p]]);
			_neighbors[p] = linked.neighbors[permutation[p]];
		}
		send.swap(_send);
		recv.swap(_recv);
		linked.neighbors.swap(_neighbors);
	};

	struct __unode__ { // unknown node
		esint offset;
		int rank;

		bool operator<(const __unode__ &other) { return rank == other.rank ? offset < other.offset : rank < other.rank; }
	};

	std::vector<__unode__> unknown;
	{ // compute size of found nodes -> if the size is lower than size of requests, we need to ask for unknown nodes
		size_t size = 0;
		for (size_t r = 0; r < recv.size(); ++r) {
			size += recv[r].front();
		}
		if (size < adjacent.size()) {
			std::vector<esint> found;
			found.reserve(size);
			unknown.reserve(adjacent.size() - size);
			for (size_t r = 0; r < recv.size(); ++r) {
				found.insert(found.end(), recv[r].begin() + 1, recv[r].begin() + 1 + recv[r].front());
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
	eslog::param("UnknownNodes", unknown.size());
	eslog::ln();

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
		auto it = std::lower_bound(merged.offsets.begin(), merged.offsets.end(), unknown[i].offset);
		if (it != merged.offsets.end() && *it == unknown[i].offset) {
			size_t n = it - merged.offsets.begin();
			found.push_back(__fnode__{unknown[i], info::mpi::rank, merged.coordinates[n]});
			if (linked.neighbors.back() != unknown[i].rank) {
				linked.neighbors.push_back(unknown[i].rank);
				send.push_back({ 0 });
				recv.push_back({ 0 });
			}
			++send.back().front();
			send.back().push_back(unknown[i].offset);
		}
	}
	utils::clearVector(unknown);

	if (found.size()) {
		permute();
	}

	if (!Communication::allGatherUnknownSize(found)) {
		eslog::internalFailure("exchange found nodes.\n");
	}
	std::sort(found.begin(), found.end());

	if (found.size()) { // insert surely found nodes
		auto begin = std::lower_bound(found.begin(), found.end(), __unode__{0, info::mpi::rank}), _begin = begin;
		auto end = begin;
		auto rank = linked.neighbors.begin();
		auto rankEnd = linked.neighbors.end();
		auto ranks = linked.neighbors.size();
		while (begin != found.end() && begin->rank == info::mpi::rank) {
			while (end != found.end() && end->rank == info::mpi::rank && begin->holder == end->holder) { ++end; }
			while (rank != rankEnd && *rank < begin->holder) { ++rank; }
			std::vector<esint> *_recv;
			if (rank != rankEnd && *rank == begin->holder) {
				_recv = &recv[rank - linked.neighbors.begin()];
				_recv->front() = (esint)(end - begin);
			} else {
				linked.neighbors.push_back(begin->holder);
				rank = linked.neighbors.begin();
				rankEnd = rank + ranks;
				send.push_back({ 0 });
				recv.push_back({ (esint)(end - begin) });
				_recv = &recv.back();
			}

			for (auto it = begin; it != end; ++it) {
				_recv->push_back(it->offset);
			}
			_recv->resize(_recv->size() + utils::reinterpret_size<esint, _Point<esfloat> >(_recv->front()));
			char *pbuffer = reinterpret_cast<char*>(_recv->data() + _recv->front() + 1);
			for (auto it = begin; it != end; ++it, pbuffer += sizeof(_Point<esfloat>)) {
				memcpy(pbuffer, &it->coordinate, sizeof(_Point<esfloat>));
			}
			begin = end;
		}
		if (_begin != end) {
			permute();
		}
	}
	eslog::checkpointln("LINKUP: UNKNOWN NODES INCLUDED");
}

// this method can be optimized by requesting nodes to closest buckets only
void linkup(MergedNodes &merged, LinkedNodes &linked, ClusteredElements &elements)
{
	eslog::startln("LINKUP: CONNECTING CLUSTERS", "LINKUP");

	// build set of required nodes (nodes that are needed by my elements)
	// nodes do not reside on this process are collected to adjacent nodes are requested by SFC neighbors
	// since SFC was used, majority of required nodes should be local
	//
	// one instance of 'adjacent' can be probably bottleneck for more than 10k MPI processes
	// hence, computation can be improved by computing adjacent to each process separately (based on SFC)

	ivector<esint> required; required.reserve(elements.enodes.size());
	for (size_t e = 0, eoffset = 0; e < elements.offsets.size(); eoffset += Element::encode(elements.etype[e++]).nodes) {
		PolyElement poly(elements.etype[e], elements.enodes.data() + eoffset);
		for (int n = 0; n < Element::encode(elements.etype[e]).nodes; ++n) {
			if (poly.isNode(n)) {
				required.push_back(elements.enodes[n + eoffset]);
			}
		}
	}
	std::vector<esint> adjacent; // nodes held by other processes
	utils::sortAndRemoveDuplicates(required);
	for (size_t offset = 0, node = 0; offset < merged.offsets.size() || node < required.size(); ++offset) {
		while (node < required.size() && (offset == merged.offsets.size() || required[node] < merged.offsets[offset])) {
			adjacent.push_back(required[node++]);
		}
		if (node < required.size() && required[node] == merged.offsets[offset]) {
			++node;
		}
	}

	eslog::checkpointln("LINKUP: ADJACENT NODES COMPUTED");
	eslog::param("AdjacentNodes", adjacent.size());
	eslog::ln();

	std::vector<std::vector<esint> > requested(linked.neighbors.size());
	if (!Communication::exchangeUnknownSize(adjacent, requested, linked.neighbors)) {
		eslog::internalFailure("request for coordinates.\n");
	}
	eslog::checkpointln("LINKUP: ADJACENT NODES EXCHANGED");

	// all nodes that were send to neighboring processes are collected in send
	// all nodes that were received from neighboring processes are collected in received
	//
	// structure:
	// send[neighbor][0] = node count N
	// send[neighbor][1..N] = node offsets
	// send[neighbor][N..2N] = node coordinates
	std::vector<std::vector<esint> > send(linked.neighbors.size()), received(linked.neighbors.size());
	for (size_t r = 0; r < requested.size(); r++) {
		send[r].push_back(0);
		for (size_t n = 0; n < requested[r].size(); n++) {
			auto nit = std::lower_bound(merged.offsets.begin(), merged.offsets.end(), requested[r][n]);
			if (nit != merged.offsets.end() && *nit == requested[r][n]) {
				send[r].push_back(nit - merged.offsets.begin());
			}
		}
		send[r].front() = send[r].size() - 1;
		send[r].resize(send[r].size() + utils::reinterpret_size<esint, _Point<esfloat> >(send[r].front()));
		char *pbuffer = reinterpret_cast<char*>(send[r].data() + send[r].front() + 1);
		for (esint n = 0; n < send[r].front(); ++n, pbuffer += sizeof(_Point<esfloat>)) {
			memcpy(pbuffer, merged.coordinates.data() + send[r][n + 1], sizeof(_Point<esfloat>));
			send[r][n + 1] = merged.offsets[send[r][n + 1]];
		}
	}
	eslog::checkpointln("LINKUP: NODES REQUESTS PROCESSED");

	if (!Communication::exchangeUnknownSize(send, received, linked.neighbors)) {
		eslog::internalFailure("return requested IDs.\n");
	}
	eslog::checkpointln("LINKUP: ADJACENT NODES RETURNED");

	// some nodes are not on SFC neighbors
	// we need to search those nodes and add update the list of neighboring processes
	exchangeUnknown(merged, linked, adjacent, send, received);
	utils::clearVector(adjacent);

	// build ranks where node resides according to requests from other processes
	//
	// there is one instance of each node in merged nodes
	// hence, all processes have to ask the holder, then holder can compose ranks information
	// warning: only the holder is aware of the node duplication
	//
	// with ranks information, we have to exchange also duplication data
	// in order to be able project all duplicated nodes to the origin
	std::vector<esint> rankDistribution(merged.offsets.size() + 1), duplicationIndex(merged.offsets.size(), -1);
	for (size_t i = 0, j = 0; i < required.size(); ++i) { // required size is proportional to clustered.noffsets
		while (j < merged.offsets.size() && merged.offsets[j] < required[i]) { ++j; }
		if (j < merged.offsets.size() && merged.offsets[j] == required[i]) {
			++rankDistribution[j];
		}
	}
	for (size_t r = 0; r < send.size(); ++r) {
		for (esint i = 1; i <= send[r].front(); ++i) { // send[r] is much smaller than clustered.noffsets
			size_t n = std::lower_bound(merged.offsets.begin(), merged.offsets.end(), send[r][i]) - merged.offsets.begin();
			++rankDistribution[n];
		}
	}
	for (size_t i = 0; i < merged.duplication.size(); ++i) {
		// real number of ranks can be lower due to the same ranks in different duplicated nodes
		rankDistribution[merged.duplication[i].origin] += rankDistribution[merged.duplication[i].duplicate];
		duplicationIndex[merged.duplication[i].origin] = i;      // always last occurrence
		duplicationIndex[merged.duplication[i].duplicate] = i; // always pointer to offset
	}
	for (size_t i = 0; i < merged.duplication.size(); ++i) {
		rankDistribution[merged.duplication[i].duplicate] = rankDistribution[merged.duplication[i].origin];
	}
	utils::sizesToOffsets(rankDistribution);
	std::vector<int, initless_allocator<int> > rankData(rankDistribution.back());

	{ // put all ranks to the rankData
		auto _rankDistribution = rankDistribution;
		size_t n = 0;
		for ( ; n < linked.neighbors.size() && linked.neighbors[n] < info::mpi::rank; ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				size_t ni = std::lower_bound(merged.offsets.begin(), merged.offsets.end(), send[n][i]) - merged.offsets.begin();
				rankData[_rankDistribution[ni]++] = linked.neighbors[n];
			}
		}
		for (size_t i = 0, j = 0; i < required.size(); ++i) { // required size is proportional to clustered.noffsets
			while (j < merged.offsets.size() && merged.offsets[j] < required[i]) { ++j; }
			if (j < merged.offsets.size() && merged.offsets[j] == required[i]) {
				rankData[_rankDistribution[j]++] = info::mpi::rank;
			}
		}
		for ( ; n < linked.neighbors.size(); ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				size_t ni = std::lower_bound(merged.offsets.begin(), merged.offsets.end(), send[n][i]) - merged.offsets.begin();
				rankData[_rankDistribution[ni]++] = linked.neighbors[n];
			}
		}

		for (size_t i = 0; i < merged.duplication.size(); ++i) {
			const DataDuplication &dd = merged.duplication[i];
			for (esint k = rankDistribution[dd.origin], j = rankDistribution[dd.duplicate]; j < _rankDistribution[dd.duplicate]; ++j) {
				while (k < _rankDistribution[dd.origin] && rankData[k] < rankData[j]) { ++k; }
				if (k == _rankDistribution[dd.origin] || rankData[k] != rankData[j]) {
					rankData[_rankDistribution[dd.origin]++] = rankData[j];
				}
			}
			std::sort(rankData.begin() + rankDistribution[dd.origin], rankData.begin() + _rankDistribution[dd.origin]);
		}
		for (size_t i = 0; i < merged.duplication.size(); ++i) {
			const DataDuplication &dd = merged.duplication[i];
			for (esint j = rankDistribution[dd.origin], k = rankDistribution[dd.duplicate]; j < _rankDistribution[dd.origin]; ++j, ++k) {
				rankData[k] = rankData[j];
			}
			_rankDistribution[dd.duplicate] = rankDistribution[dd.duplicate] + _rankDistribution[dd.origin] - rankDistribution[dd.origin];
		}

		size_t sum = 0;
		for (size_t i = 0, j = 0; i < merged.offsets.size(); ++i, sum = j) {
			for (esint k = rankDistribution[i]; k < _rankDistribution[i]; ++k) {
				rankData[j++] = rankData[k];
			}
			rankDistribution[i] = sum;
		}
		rankDistribution.back() = sum;
		rankData.resize(rankDistribution.back());
	}

	// compute list of ranks in the form:
	// ranks[neighbors]{ ranks count; duplication count; ranks; duplications; ... }
	std::vector<std::vector<esint> > ranks(linked.neighbors.size());
	{ // send rank data and duplication info
		std::vector<std::vector<esint> > sBuffer(linked.neighbors.size());

		for (size_t n = 0; n < send.size(); ++n) {
			size_t size = 0;
			for (esint i = 1; i <= send[n].front(); ++i) {
				size += 2; // ranks, duplications
				esint ni = std::lower_bound(merged.offsets.begin(), merged.offsets.end(), send[n][i]) - merged.offsets.begin();
				size += rankDistribution[ni + 1] - rankDistribution[ni];
				if (duplicationIndex[ni] != -1) {
					if (ni == merged.duplication[duplicationIndex[ni]].duplicate) {
						ni = merged.duplication[duplicationIndex[ni]].origin;
					}
					esint di = duplicationIndex[ni];
					while (di >= 0 && ni == merged.duplication[di].origin) {
						size += 2;
						--di;
					}
				}
			}
			sBuffer[n].reserve(size);
		}
		for (size_t n = 0; n < send.size(); ++n) {
			for (esint i = 1; i <= send[n].front(); ++i) {
				esint ni = std::lower_bound(merged.offsets.begin(), merged.offsets.end(), send[n][i]) - merged.offsets.begin();
				size_t prevsize = sBuffer[n].size();
				sBuffer[n].push_back(rankDistribution[ni + 1] - rankDistribution[ni]);
				sBuffer[n].push_back(0);
				for (esint r = rankDistribution[ni]; r < rankDistribution[ni + 1]; ++r) {
					sBuffer[n].push_back(rankData[r]);
				}
				if (duplicationIndex[ni] != -1) {
					if (ni == merged.duplication[duplicationIndex[ni]].duplicate) {
						ni = merged.duplication[duplicationIndex[ni]].origin;
					}
					esint di = duplicationIndex[ni];
					while (di >= 0 && ni == merged.duplication[di].origin) {
						++sBuffer[n][prevsize + 1];
						sBuffer[n].push_back(merged.offsets[merged.duplication[di].origin]);
						sBuffer[n].push_back(merged.offsets[merged.duplication[di].duplicate]);
						--di;
					}
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sBuffer, ranks, linked.neighbors)) {
			eslog::internalFailure("cannot exchange clustered nodes ranks.\n");
		}
	}

	{ // build nodes
		// required nodes are not projected
		// duplication info is not local anymore -> local offsets are changed to global offsets

		struct __node__ {
			esint offset, nranks, *ranks;
			const _Point<esfloat>* coordinate;

			bool operator<(const __node__ &other) const { return other.offset < offset; }
		};

		std::priority_queue<__node__> projection;
		std::unordered_set<int> neighbors;
		linked.offsets.reserve(required.size());
		linked.coordinates.reserve(required.size());
		linked.rankDistribution.reserve(required.size() + 1);
		linked.rankDistribution.push_back(0);
		for (auto dup = merged.duplication.begin(); dup != merged.duplication.end(); ++dup) {
			linked.duplication.push_back({ merged.offsets[dup->origin], merged.offsets[dup->duplicate] }); // it can include unnecessary nodes
		}

		auto push = [&] (const esint &offset, const _Point<esfloat> *coo, const esint &nranks, const int *ranks) {
			linked.offsets.push_back(offset);
			linked.coordinates.push_back(*coo);
			for (esint n = 0; n < nranks; ++n) {
				linked.rankData.push_back(ranks[n]);
				if (linked.rankData.back() != info::mpi::rank) {
					neighbors.insert(linked.rankData.back());
				}
			}
			linked.rankDistribution.push_back(linked.rankData.size());
		};

		std::vector<esint> recvOffset(linked.neighbors.size()), rankOffset(linked.neighbors.size());
		for (size_t i = 0, j = 0; i < required.size(); ++i) {
			while (projection.size() && projection.top().offset <= required[i]) { // it works since duplicate < origin
				while (projection.size() && projection.top().offset <= required[i]) {
					if (projection.top().offset < required[i]) {
						push(projection.top().offset, projection.top().coordinate, projection.top().nranks, projection.top().ranks);
						while (projection.size() && linked.offsets.back() == projection.top().offset) { projection.pop(); }
					} else {
						while (projection.size() && required[i] == projection.top().offset) { projection.pop(); }
					}
				}
			}
			while (j < merged.offsets.size() && merged.offsets[j] < required[i]) { ++j; }
			if (j < merged.offsets.size() && merged.offsets[j] == required[i]) { // local node
				if (duplicationIndex[j] != -1 && merged.offsets[j] == linked.duplication[duplicationIndex[j]].duplicate) {
					projection.push({
						linked.duplication[duplicationIndex[j]].origin,
						rankDistribution[j + 1] - rankDistribution[j],
						rankData.data() + rankDistribution[j],
						&merged.coordinates[j]});
				} else {
					push(merged.offsets[j], &merged.coordinates[j], rankDistribution[j + 1] - rankDistribution[j], rankData.data() + rankDistribution[j]);
				}
			} else { // adjacent node
				for (size_t r = 0; r < recvOffset.size(); ++r) {
					if (recvOffset[r] < received[r].front() && required[i] == received[r][recvOffset[r] + 1]) {
						const _Point<esfloat> *coo = reinterpret_cast<_Point<esfloat>*>(received[r].data() + received[r].front() + 1) + recvOffset[r]++;
						esint nranks = ranks[r][rankOffset[r]++];
						esint nduplications = ranks[r][rankOffset[r]++];
						esint project = -1;
						for (esint n = 0; n < nduplications; ++n) {
							linked.duplication.push_back({ ranks[r][rankOffset[r] + nranks + 2 * n], ranks[r][rankOffset[r] + nranks + 2 * n + 1] });
							if (linked.duplication.back().duplicate == required[i]) {
								project = linked.duplication.back().origin;
							}
						}
						if (project != -1) {
							projection.push({ project, nranks, ranks[r].data() + rankOffset[r], coo});
						} else {
							push(required[i], coo, nranks, ranks[r].data() + rankOffset[r]);
						}
						rankOffset[r] += nranks + 2 * nduplications;
					}
				}
			}
		}
		while (projection.size()) {
			push(projection.top().offset, projection.top().coordinate, projection.top().nranks, projection.top().ranks);
			while (projection.size() && linked.offsets.back() == projection.top().offset) { projection.pop(); }
		}
		utils::sortAndRemoveDuplicates(linked.duplication);

		linked.neighbors.clear();
		for (auto n = neighbors.begin(); n != neighbors.end(); ++n) {
			linked.neighbors.push_back(*n);
		}
		std::sort(linked.neighbors.begin(), linked.neighbors.end());
	}

	eslog::endln("LINKUP: LINKED UP");
}

}
}
