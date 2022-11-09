
#include "builder.utils.h"

#include "basis/containers/tarray.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/structures/kdtree.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "wrappers/mpi/communication.h"

#include <array>
#include <functional>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace espreso {
namespace builder {

static void treeSearch(std::vector<_Point<esfloat>, initless_allocator<_Point<esfloat> > > &coordinates, std::function<void(esint origin, esint target)> merge)
{
	KDTree<esfloat> tree(coordinates.data(), coordinates.data() + coordinates.size());

	profiler::synccheckpoint("kdtree_build");
	profiler::param("kdtree_levels", tree.levels);

	double eps = info::ecf->input.duplication_tolerance;
	std::vector<esint> duplicate(tree.permutation.size(), -1);

	for (esint i = std::exp2(tree.levels), first = i; i < std::exp2(tree.levels + 1); ++i) {
		esint begin = tree.begin(i);
		esint end = tree.end(i);
		if (begin == end) {
			continue;
		}

		_Point<esfloat> min;
		tree.boxMin(i, min);

		auto check = [&] (esint p, esint begin, esint end) {
			for (auto pp = tree.permutation.cbegin() + begin; pp != tree.permutation.cbegin() + end; ++pp) {
				if (
						coordinates[tree.permutation[p]].x <= coordinates[*pp].x + eps && coordinates[*pp].x - eps <= coordinates[tree.permutation[p]].x &&
						coordinates[tree.permutation[p]].y <= coordinates[*pp].y + eps && coordinates[*pp].y - eps <= coordinates[tree.permutation[p]].y &&
						coordinates[tree.permutation[p]].z <= coordinates[*pp].z + eps && coordinates[*pp].z - eps <= coordinates[tree.permutation[p]].z) {

					if (duplicate[pp - tree.permutation.cbegin()] >= 0) {
						merge(tree.permutation[p], tree.permutation[duplicate[pp - tree.permutation.cbegin()]]);
						duplicate[p] = duplicate[pp - tree.permutation.cbegin()];
					} else {
						merge(tree.permutation[p], *pp);
						duplicate[p] = pp - tree.permutation.cbegin();
					}
					break;
				}
			}
		};

		std::function<void(size_t, size_t, esint)> traverse = [&] (size_t node, size_t max, esint p) {
			if (coordinates[tree.permutation[p]][tree.splitters[node].d] <= tree.splitters[node].value + eps) {
				if (2 * node < tree.splitters.size()) {
					traverse(2 * node, max, p);
				} else {
					if (2 * node < max) {
						check(p, tree.begin(2 * node), tree.end(2 * node));
					}
				}
			}
			if (tree.splitters[node].value - eps <= coordinates[tree.permutation[p]][tree.splitters[node].d]) {
				if (2 * node < tree.splitters.size()) {
					traverse(2 * node + 1, max, p);
				} else {
					if (2 * node + 1 < max) {
						check(p, tree.begin(2 * node + 1), tree.end(2 * node + 1));
					}
				}
			}
		};

		if (tree.splitters.size() > 1) {
			for (auto p = tree.permutation.cbegin() + begin; p != tree.permutation.cbegin() + end; ++p) {
				if ((coordinates[*p].x <= min.x + eps) || (coordinates[*p].y <= min.y + eps) || (coordinates[*p].z <= min.z + eps)) {
					traverse(1, std::max(first, i), p - tree.permutation.cbegin());
				}
			}
		}

		for (auto left = tree.permutation.cbegin() + begin, right = left + 1; right != tree.permutation.cbegin() + end; ++right) {
			if (duplicate[right - tree.permutation.cbegin()] != -1) {
				continue;
			}
			while (left != right && (coordinates[*left].x + eps < coordinates[*right].x || duplicate[left - tree.permutation.cbegin()] != -1)) {
				++left;
			}
			for (auto mid = left; mid != right;) {
				while (mid != right && coordinates[*mid].y + eps <  coordinates[*right].y) {
					++mid;
				}
				if (mid != right && coordinates[*mid].y - eps <= coordinates[*right].y) {
					if (duplicate[mid - tree.permutation.cbegin()] == -1 && coordinates[*right].z <= coordinates[*mid].z + eps && coordinates[*mid].z - eps <= coordinates[*right].z) {
						if (duplicate[mid - tree.permutation.cbegin()] >= 0) {
							merge(*right, tree.permutation[duplicate[mid - tree.permutation.cbegin()]]);
							duplicate[right - tree.permutation.cbegin()] = duplicate[mid - tree.permutation.cbegin()];
						} else {
							merge(*right, *mid);
							duplicate[right - tree.permutation.cbegin()] = mid - tree.permutation.cbegin();
						}
					}
					++mid;
				}
				while (mid != right && coordinates[*right].y + eps < coordinates[*mid].y) {
					++mid;
				}
			}
		}
	}
}

void searchDuplicatedNodes(ClusteredNodes &clustered, MergedNodes &merged)
{
	merged.coordinates.swap(clustered.coordinates);
	merged.offsets.swap(clustered.offsets);

	treeSearch(merged.coordinates, [&] (esint origin, esint target) {
		merged.duplication.push_back({ target, origin });
	});
	std::sort(merged.duplication.begin(), merged.duplication.end());

	// set origin to the last occurrence -> it simplify linking step
	for (auto begin = merged.duplication.begin(), end = begin; begin != merged.duplication.end(); begin = end) {
		while (end != merged.duplication.end() && begin->origin == end->origin) { ++end; }
		esint max = std::max((end - 1)->origin, (end - 1)->duplicate);
		for (auto it = begin; it != end; ++it) {
			if (it->duplicate == max) {
				std::swap(it->duplicate, it->origin);
			} else {
				it->origin = max;
			}
		}
	}
	std::sort(merged.duplication.begin(), merged.duplication.end());

	size_t duplicates = merged.duplication.size();
	Communication::allReduce(&duplicates, nullptr, 1, MPITools::getType(duplicates).mpitype, MPI_SUM);
	eslog::info(" == DUPLICATED NODES %70lu == \n", duplicates);
}

void mergeDuplicatedNodes(MergedNodes &merged)
{
	ivector<esint> usedNodes(merged.coordinates.size());
	std::iota(usedNodes.begin(), usedNodes.end(), 0);
	for (size_t n = 0; n < merged.duplication.size(); ++n) {
		usedNodes[merged.duplication[n].duplicate] = merged.duplication[n].origin;
	}

	for (size_t n = 0, offset = 0; n < usedNodes.size(); ++n) {
		if (usedNodes[n] == (esint)n) {
			usedNodes[n] = offset++;
		}
	}

	size_t offset = 0;
	for (size_t n = 0; n < usedNodes.size(); ++n) {
		if (usedNodes[n] <= (esint)n) {
			merged.offsets[offset] = merged.offsets[n];
			merged.coordinates[offset] = merged.coordinates[n];
			++offset;
		} else {
			usedNodes[n] = usedNodes[usedNodes[n]];
		}
	}
	merged.offsets.resize(offset);
	merged.coordinates.resize(offset);
}

void exchangeSFCBoundaryNodes(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, const std::vector<int> &sfcNeighbors, ClusteredNodes &clustered)
{
	std::vector<esint> cdistribution = tarray<esint>::distribute(info::env::OMP_NUM_THREADS, clustered.offsets.size());

	double eps = info::ecf->input.duplication_tolerance;
	eps += info::ecf->input.duplication_tolerance * 1e-6;

	struct __node__ {
		esint offset;
		int target;

		bool operator<(const __node__ &other) { return target == other.target ? offset < other.offset : target < other.target; }
	};

	std::vector<std::vector<__node__> > toneighs(info::env::OMP_NUM_THREADS);
	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<__node__> tneighs;
		std::vector<int> neighs;
		neighs.reserve(27);
		for (esint n = cdistribution[t]; n < cdistribution[t + 1]; ++n) {
			int hit[6] = { 0, 0, 0, 0, 0, 0};
			for (size_t d = 0; d < sfc.dimension; ++d) {
				size_t origin = sfc.getBucket(clustered.coordinates[n][d], d);
				if (sfc.getBucket(clustered.coordinates[n][d] - eps, d) < origin) {
					hit[2 * d + 0] = -1;
				}
				if (sfc.getBucket(clustered.coordinates[n][d] + eps, d) > origin) {
					hit[2 * d + 1] = +1;
				}
			}

			neighs.clear();
			for (int x = hit[0]; x <= hit[1]; ++x) {
				for (int y = hit[2]; y <= hit[3]; ++y) {
					for (int z = hit[4]; z <= hit[5]; ++z) {
						if (x || y || z) {
							size_t b = sfc.getBucket(clustered.coordinates[n] + Point(x * eps, y * eps, z * eps));
							int nn = std::lower_bound(splitters.begin(), splitters.end(), b + 1) - splitters.begin() - 1;
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
				tneighs.push_back(__node__{n, *nn});
			}
		}
		toneighs[t].swap(tneighs);
	}

	{ // merge and sort
		size_t size = 0;
		for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
			size += toneighs[t].size();
		}
		toneighs.reserve(size);

		for (int t = 1; t < info::env::OMP_NUM_THREADS; t++) {
			toneighs[0].insert(toneighs[0].end(), toneighs[t].begin(), toneighs[t].end());
			utils::clearVector(toneighs[t]);
		}
		std::sort(toneighs[0].begin(), toneighs[0].end());
	}

	std::vector<esint> offsets;
	ivector<_Point<esfloat> > coordinates;
	std::vector<std::vector<esint> > sBuffer(sfcNeighbors.size()), rBuffer(sfcNeighbors.size());
	{ // build send buffer
		size_t ni = 0;
		for (auto begin = toneighs[0].begin(), end = begin; end != toneighs[0].end(); begin = end) {
			while (sfcNeighbors[ni] < begin->target) { ++ni; }
			while (end != toneighs[0].end() && begin->target == end->target) { ++end; }

			for (auto n = begin; n != end; ++n) {
				offsets.push_back(n->offset);
				coordinates.push_back(clustered.coordinates[n->offset]);
			}
			esint nodes = end - begin;
			sBuffer[ni].resize(1 + nodes + utils::reinterpret_size<esint, _Point<esfloat> >(nodes));
			sBuffer[ni][0] = nodes;
			memcpy(sBuffer[ni].data() + 1, offsets.data() + offsets.size() - nodes, nodes * sizeof(esint));
			memcpy(sBuffer[ni].data() + 1 + nodes, reinterpret_cast<esint*>(coordinates.data() + coordinates.size() - nodes), nodes * sizeof(_Point<esfloat>));
		}
	}
	esint mynodes = offsets.size();
	utils::clearVector(toneighs[0]);

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, sfcNeighbors)) {
		eslog::internalFailure("cannot exchange potentially duplicated nodes.\n");
	}
	utils::clearVector(sBuffer);

	{ // preallocate data
		size_t size = offsets.size();
		for (size_t n = 0; n < rBuffer.size(); ++n) {
			if (rBuffer[n].size()) {
				size += rBuffer[n].front();
			}
		}
		offsets.reserve(size);
		coordinates.reserve(size);
	}

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		if (rBuffer[n].size()) {
			esint nodes = rBuffer[n].front();
			offsets.insert(offsets.end(), rBuffer[n].begin() + 1, rBuffer[n].begin() + 1 + nodes);
			coordinates.insert(coordinates.end(), reinterpret_cast<_Point<esfloat>*>(rBuffer[n].data() + 1 + nodes), reinterpret_cast<_Point<esfloat>*>(rBuffer[n].data() + 1 + nodes) + nodes);
		}
	}

	std::vector<size_t> toinsert, toerase;
	treeSearch(coordinates, [&] (esint origin, esint target) {
		if (origin < mynodes && mynodes <= target) { // I have duplicate
			toerase.push_back(offsets[origin]);
		}
		if (target < mynodes && mynodes <= origin) { // rank n has duplicate and I have target
			toinsert.push_back(origin);
		}
		// rest of nodes are merged later
	});

	if (toerase.size()) {
		utils::sortAndRemoveDuplicates(toerase);
	}
	if (toinsert.size()) {
		utils::sortAndRemoveDuplicates(toinsert);
	}

	size_t size = clustered.offsets.size();
	if (toerase.size()) {
		for (size_t n = toerase.front(), last = toerase.front(), e = 0; n < size; ++n) {
			if (toerase.size() <= e || n != toerase[e]) {
				clustered.offsets[last] = clustered.offsets[n];
				clustered.coordinates[last] = clustered.coordinates[n];
				++last;
			} else {
				++e;
			}
		}
		size -= toerase.size();
		clustered.coordinates.resize(size);
		clustered.offsets.resize(size);
	}

	if (toinsert.size()) {
		clustered.coordinates.resize(size + toinsert.size());
		clustered.offsets.resize(size + toinsert.size());
		for (size_t i = size, j = 0; j < toinsert.size(); ++i, ++j) {
			clustered.offsets[i] = offsets[toinsert[j]];
			clustered.coordinates[i] = coordinates[toinsert[j]];
		}
	}
}

void globalToLocal(ClusteredElements &clustered, MergedElements &merged, LinkedNodes &nodes)
{
	// duplicate boundary is adept to duplication only
	std::unordered_map<esint, esint> g2l;

	g2l.reserve(nodes.offsets.size() + nodes.duplication.size());
	for (size_t n = 0; n < nodes.offsets.size(); ++n) {
		g2l[nodes.offsets[n]] = n;
	}
	for (auto dup = nodes.duplication.begin(); dup != nodes.duplication.end(); ++dup) {
		g2l[dup->duplicate] = g2l[dup->origin];
	}

	{ // build elements
		swap(merged, clustered);
		merged.offsets.swap(clustered.offsets);
	}

	{ // build nodes
		std::unordered_map<esint, esint> dmap;
		for (size_t n = 0; n < nodes.duplication.size(); ++n) {
			dmap[nodes.duplication[n].duplicate] = g2l[nodes.duplication[n].origin];
		}
		std::unordered_map<esint, esint> rmap;
		for (size_t n = 0; n < nodes.neighbors.size(); ++n) {
			rmap[nodes.neighbors[n]] = n;
		}
		std::vector<int> neighbors(nodes.neighbors.size());

		std::vector<esint> usedNode(nodes.offsets.size() + 1);
		for (size_t e = 0, eoffset = 0; e < merged.etype.size(); eoffset += Element::encode(merged.etype[e++]).nodes) {
			PolyElement poly(merged.etype[e], merged.enodes.data() + eoffset);
			for (int n = 0; n < Element::encode(merged.etype[e]).nodes; ++n) {
				if (poly.isNode(n)) {
					usedNode[g2l[merged.enodes[n + eoffset]]] = 1;
				}
			}
		}
		utils::sizesToOffsets(usedNode);
		for (size_t e = 0, eoffset = 0; e < merged.etype.size(); eoffset += Element::encode(merged.etype[e++]).nodes) {
			PolyElement poly(merged.etype[e], merged.enodes.data() + eoffset);
			for (int n = 0; n < Element::encode(merged.etype[e]).nodes; ++n) {
				if (poly.isNode(n)) {
					merged.enodes[n + eoffset] = usedNode[g2l[merged.enodes[n + eoffset]]];
				}
			}
		}

		// we need to inform about erased nodes
		std::vector<std::vector<esint> > sRemoved(nodes.neighbors.size()), rRemoved(nodes.neighbors.size());
		for (size_t n = 0; n < nodes.offsets.size(); ++n) {
			if (usedNode[n] == usedNode[n + 1] && usedNode[g2l[nodes.offsets[n]]] == usedNode[g2l[nodes.offsets[n]] + 1]) {
				for (esint r = nodes.ranks.distribution[n]; r < nodes.ranks.distribution[n + 1]; ++r) {
					if (nodes.ranks.data[r] != info::mpi::rank) {
						sRemoved[rmap[nodes.ranks.data[r]]].push_back(nodes.offsets[n]);
					}
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sRemoved, rRemoved, nodes.neighbors)) {
			eslog::internalFailure("exchange removed nodes.\n");
		}

		std::vector<std::pair<esint, esint> > removed;
		for (size_t r = 0; r < rRemoved.size(); ++r) {
			for (size_t i = 0; i < rRemoved[r].size(); ++i) {
				auto dup = dmap.find(rRemoved[r][i]);
				if (dup != dmap.end()) {
					removed.push_back(std::pair<esint, esint>(dup->second, nodes.neighbors[r]));
				} else {
					removed.push_back(std::pair<esint, esint>(rRemoved[r][i], nodes.neighbors[r]));
				}
			}
		}
		std::sort(removed.begin(), removed.end());

		size_t offset = 0, ranks = 0, rindex = 0;
		for (size_t n = 0; n < nodes.offsets.size(); ++n) {
			if (usedNode[n] < usedNode[n + 1]) {
				while (rindex < removed.size() && removed[rindex].first < nodes.offsets[n]) { ++rindex; }
				nodes.offsets[offset] = nodes.offsets[n];
				nodes.coordinates[offset] = nodes.coordinates[n];
				size_t rdist = ranks;
				for (esint r = nodes.ranks.distribution[n]; r < nodes.ranks.distribution[n + 1]; ++r) {
					if (rindex == removed.size() || removed[rindex].first != nodes.offsets[n] || removed[rindex].second != nodes.ranks.data[r]) {
						nodes.ranks.data[ranks++] = nodes.ranks.data[r];
						if (nodes.ranks.data[r] != info::mpi::rank) {
							++neighbors[rmap[nodes.ranks.data[r]]];
						}
					}
					if (rindex < removed.size() && removed[rindex].first == nodes.offsets[n] && removed[rindex].second == nodes.ranks.data[r]) {
						++rindex;
					}
				}
				nodes.ranks.distribution[offset++] = rdist;
			}
		}
		nodes.offsets.resize(offset);
		nodes.coordinates.resize(offset);
		nodes.ranks.distribution.resize(offset + 1);
		nodes.ranks.distribution.back() = ranks;
		nodes.ranks.data.resize(ranks);
		std::vector<int> neighs; neighs.swap(nodes.neighbors);
		for (size_t n = 0; n < neighbors.size(); ++n) {
			if (neighbors[n]) {
				nodes.neighbors.push_back(neighs[n]);
			}
		}
	}

	eslog::info(" == DUPLICATED ELEMENTS %67lu == \n", 0);
}

void mergeDuplicatedElements(ClusteredElements &clustered, MergedElements &merged, LinkedNodes &nodes, int dimension)
{
	// duplicate boundary is adept to duplication only
	std::unordered_set<esint> duplicateBoundary, duplicateElement;
	std::unordered_map<esint, esint> g2l;

	g2l.reserve(nodes.offsets.size() + nodes.duplication.size());
	for (size_t n = 0; n < nodes.offsets.size(); ++n) {
		g2l[nodes.offsets[n]] = n;
	}
	for (auto dup = nodes.duplication.begin(); dup != nodes.duplication.end(); ++dup) {
		g2l[dup->duplicate] = g2l[dup->origin];
	}

	// eoffset, etype, enodes; .... + project all enodes
	std::vector<std::vector<esint> > sBuffer(nodes.neighbors.size()), rBuffer(nodes.neighbors.size());
	if (nodes.neighbors.size()) { // duplicated elements have to have all nodes held by other rank
		std::vector<int> neighbors;
		for (size_t e = 0, eoffset = 0; e < clustered.offsets.size(); eoffset += Element::encode(clustered.etype[e++]).nodes) {
			neighbors.assign(nodes.neighbors.begin(), nodes.neighbors.end());
			PolyElement poly(clustered.etype[e], clustered.enodes.data() + eoffset);
			for (int n = 0; n < Element::encode(clustered.etype[e]).nodes; ++n) {
				if (poly.isNode(n)) {
					esint local = g2l[clustered.enodes[n + eoffset]];
					size_t intersection = 0, current = 0;
					for (int r = nodes.ranks.distribution[local]; r < nodes.ranks.distribution[local + 1]; ++r) {
						while (current < neighbors.size() && neighbors[current] < nodes.ranks.data[r]) { ++current; }
						if (current < neighbors.size() && neighbors[current] == nodes.ranks.data[r]) {
							neighbors[intersection++] = neighbors[current++];
						}
					}
					neighbors.resize(intersection);
				}
			}
			esint roffset = 0;
			for (auto r = neighbors.begin(); r != neighbors.end(); ++r) {
				while (nodes.neighbors[roffset] < *r) { ++roffset; }
				sBuffer[roffset].push_back(clustered.offsets[e]);
				sBuffer[roffset].push_back((esint)clustered.etype[e]);
				for (size_t n = eoffset; n < eoffset + Element::encode(clustered.etype[e]).nodes; ++n) {
					sBuffer[roffset].push_back(clustered.enodes[n]);
				}
			}
			if (neighbors.size() && Mesh::element(clustered.etype[e]).dimension < dimension) {
				duplicateBoundary.insert(e);
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, nodes.neighbors)) {
		eslog::internalFailure("cannot exchange duplicated elements.\n");
	}
	utils::clearVector(sBuffer);

	auto getMinimal = [&] (const Element::CODE &code, esint *enodes) {
		PolyElement poly(code, enodes);
		int n = 0;
		while (!poly.isNode(n)) ++n;
		esint min = g2l[enodes[n]];
		for (++n; n < Element::encode(code).nodes; ++n) {
			if (poly.isNode(n)) {
				esint current = g2l[enodes[n]];
				if (nodes.coordinates[current].x < nodes.coordinates[min].x) {
					min = current;
				} else if (nodes.coordinates[current].x == nodes.coordinates[min].x && current < min) {
					min = current;
				}
			}
		}
		return min;
	};

	std::vector<esint> mapDist(nodes.offsets.size() + 1);
	ivector<esint> mapData;

	for (size_t e = 0, eoffset = 0; e < clustered.offsets.size(); eoffset += Element::encode(clustered.etype[e++]).nodes) {
		if (Mesh::element(clustered.etype[e]).dimension == dimension) {
			++mapDist[getMinimal(clustered.etype[e], clustered.enodes.data() + eoffset)];
		} else {
			if (duplicateBoundary.count(e)) { // we have to search a parent for each sent elements
				PolyElement poly(clustered.etype[e], clustered.enodes.data() + eoffset);
				for (int n = 0; n < Element::encode(clustered.etype[e]).nodes; ++n) {
					if (poly.isNode(n)) {
						++mapDist[g2l[clustered.enodes[n + eoffset]]];
					}
				}
			}
		}
	}

	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Element::encode((Element::CODE)rBuffer[r][i + 1]).nodes) {
			if (Mesh::element(rBuffer[r][i + 1]).dimension == dimension) {
				++mapDist[getMinimal((Element::CODE)rBuffer[r][i + 1], rBuffer[r].data() + i + 2)];
			} else {
				PolyElement poly((Element::CODE)rBuffer[r][i + 1], rBuffer[r].data() + i + 2);
				for (int n = 0; n < Element::encode((Element::CODE)rBuffer[r][i + 1]).nodes; ++n) {
					if (poly.isNode(n)) {
						if (g2l.count(rBuffer[r][i + 2 + n]) == 0) {
							eslog::error("unknown node was received.\n");
						}
						++mapDist[g2l[rBuffer[r][i + 2 + n]]];
					}
				}
			}
		}
	}
	utils::sizesToOffsets(mapDist);
	mapData.resize(mapDist.back(), -1);

	std::vector<std::pair<esint, esint> > recvMap;
	std::vector<esint> _mapDist = mapDist;

	for (size_t e = 0, eoffset = 0; e < clustered.offsets.size(); eoffset += Element::encode(clustered.etype[e++]).nodes) {
		if (Mesh::element(clustered.etype[e]).dimension == dimension) {
			mapData[_mapDist[getMinimal(clustered.etype[e], clustered.enodes.data() + eoffset)]++] = e;
		} else {
			if (duplicateBoundary.count(e)) {
				PolyElement poly(clustered.etype[e], clustered.enodes.data() + eoffset);
				for (int n = 0; n < Element::encode(clustered.etype[e]).nodes; ++n) {
					if (poly.isNode(n)) {
						mapData[_mapDist[g2l[clustered.enodes[n + eoffset]]]++] = e;
					}
				}
			}
		}
	}
	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Element::encode((Element::CODE)rBuffer[r][i + 1]).nodes) {
			if (Mesh::element(rBuffer[r][i + 1]).dimension == dimension) {
				mapData[_mapDist[getMinimal((Element::CODE)rBuffer[r][i + 1], rBuffer[r].data() + i + 2)]++] = clustered.offsets.size() + recvMap.size();
			} else {
				PolyElement poly((Element::CODE)rBuffer[r][i + 1], rBuffer[r].data() + i + 2);
				for (int n = 0; n < Element::encode((Element::CODE)rBuffer[r][i + 1]).nodes; ++n) {
					if (poly.isNode(n)) {
						mapData[_mapDist[g2l[rBuffer[r][i + 2 + n]]]++] = clustered.offsets.size() + recvMap.size();
					}
				}
			}
			recvMap.push_back({ (esint)r, (esint)i });
		}
	}
	utils::clearVector(_mapDist);

	ivector<esint> shared, tmpShared;
	auto merge = [&] (const Element::CODE &code, const esint *enodes) {
		PolyElement poly(code, enodes);
		int n = 0;
		while (!poly.isNode(n)) ++n;
		esint min = g2l[enodes[n]];
		shared.assign(mapData.begin() + mapDist[min], mapData.begin() + mapDist[min + 1]);
		for (++n; n < Element::encode(code).nodes; ++n) {
			if (poly.isNode(n)) {
				esint current = g2l[enodes[n]];
				tmpShared.swap(shared);
				shared.resize(tmpShared.size() + mapDist[current + 1] - mapDist[current]);
				std::merge(mapData.begin() + mapDist[current], mapData.begin() + mapDist[current + 1], tmpShared.begin(), tmpShared.end(), shared.begin());
				if (nodes.coordinates[current].x < nodes.coordinates[min].x) {
					min = current;
				} else if (nodes.coordinates[current].x == nodes.coordinates[min].x && current < min) {
					min = current;
				}
			}
		}
		return min;
	};

	struct __parent__ {
		esint offset, parent, rank;
		bool operator<(const __parent__ &other) const { return offset == other.offset ? parent < other.parent : offset < other.offset; }
	};

	std::vector<__parent__> parents;
	std::vector<std::vector<__parent__> > rParents(nodes.neighbors.size());
	std::vector<esint> _checkBuffer;
	_checkBuffer.reserve(2 * 20);
	auto areSame = [&] (const Element::CODE &code1, esint *nodes1, const Element::CODE &code2, esint *nodes2) {
		if (code1 == code2) {
			_checkBuffer.clear();
			PolyElement poly1(code1, nodes1), poly2(code2, nodes2);
			for (esint n = 0; n < Element::encode(code1).nodes; ++n) { if (poly1.isNode(n)) _checkBuffer.push_back(g2l[nodes1[n]]); }
			for (esint n = 0; n < Element::encode(code2).nodes; ++n) { if (poly2.isNode(n)) _checkBuffer.push_back(g2l[nodes2[n]]); }
			std::sort(_checkBuffer.begin(), _checkBuffer.begin() + poly1.nodes);
			std::sort(_checkBuffer.begin() + poly1.nodes, _checkBuffer.end());
			return memcmp(_checkBuffer.data(), _checkBuffer.data() + poly1.nodes, sizeof(esint) * poly1.nodes) == 0;
		}
		return false;
	};

	std::vector<esint> edist = {0};
	for (size_t e = 0; e < clustered.offsets.size(); ++e) {
		edist.push_back(edist.back() + Element::encode(clustered.etype[e]).nodes);
	}

	for (size_t e1 = 0; e1 < clustered.offsets.size(); ++e1) {
		if (Mesh::element(clustered.etype[e1]).dimension == dimension && duplicateElement.count(e1) == 0) {
			esint min = merge(clustered.etype[e1], clustered.enodes.data() + edist[e1]);
			for (esint n = mapDist[min]; n < mapDist[min + 1]; ++n) {
				size_t e2 = mapData[n];
				if (e1 < e2) { // if e2 is lower, the pair was already processed
					if (e2 < clustered.offsets.size()) { // local
						if (Mesh::element(clustered.etype[e2]).dimension == dimension && areSame(clustered.etype[e1], clustered.enodes.data() + edist[e1], clustered.etype[e2], clustered.enodes.data() + edist[e2])) {
							merged.duplication.push_back({ clustered.offsets[e1], clustered.offsets[e2] });
							duplicateElement.insert(e2);
						}
					} else { // from neighbors
						const auto &r = recvMap[e2 - clustered.offsets.size()];
						if (Mesh::element(rBuffer[r.first][r.second + 1]).dimension == dimension) {
							if (areSame(clustered.etype[e1], clustered.enodes.data() + edist[e1], (Element::CODE)rBuffer[r.first][r.second + 1], rBuffer[r.first].data() + r.second + 2)) {
								if (clustered.offsets[e1] < rBuffer[r.first][r.second]) {
									merged.duplication.push_back({ clustered.offsets[e1], rBuffer[r.first][r.second] });
									// element is removed on the neighboring process
								} else {
									merged.duplication.push_back({ rBuffer[r.first][r.second], clustered.offsets[e1] });
									duplicateElement.insert(e1);
								}
							}
						}
					}
				}
			}

			for (auto begin = shared.begin(), end = begin; end != shared.end(); begin = end) {
				while (end != shared.end() && *begin == *end) { ++end; }
				if (end - begin > 1) { // at least Line2
					size_t e2 = *begin;
					if (e2 < clustered.offsets.size()) {
						if (PolyElement(clustered.etype[e2], clustered.enodes.data() + edist[e2]).nodes == end - begin) {
							parents.push_back({ clustered.offsets[e2], clustered.offsets[e1], info::mpi::rank });
						}
					} else {
						auto &roffset = recvMap[e2 - clustered.offsets.size()];
						if (PolyElement((Element::CODE)rBuffer[roffset.first][roffset.second + 1], rBuffer[roffset.first].data() + roffset.second + 2).nodes == end - begin) {
							parents.push_back({ rBuffer[roffset.first][roffset.second], clustered.offsets[e1], info::mpi::rank });
						}
					}
				}
			}
		}
	}

	if (!Communication::exchangeUnknownSize(parents, rParents, nodes.neighbors)) {
		eslog::internalFailure("exchange parent elements.\n");
	}
	for (size_t r = 0; r < rParents.size(); ++r) {
		parents.insert(parents.end(), rParents[r].begin(), rParents[r].end());
	}
	std::sort(parents.begin(), parents.end());
	std::sort(merged.duplication.begin(), merged.duplication.end());

	// linked.elements should have enough capacity to insert some other elements (see clusterization)
	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Element::encode((Element::CODE)rBuffer[r][i + 1]).nodes) {
			if (Mesh::element(rBuffer[r][i + 1]).dimension < dimension) { // only boundary elements can be inserted
				auto p = std::lower_bound(parents.begin(), parents.end(), __parent__{ rBuffer[r][i], 0, 0 });
				if (p != parents.end() && p->offset == rBuffer[r][i] && p->rank == info::mpi::rank) {
					clustered.offsets.push_back(rBuffer[r][i]);
					clustered.etype.push_back((Element::CODE)rBuffer[r][i + 1]);
					clustered.enodes.insert(clustered.enodes.end(), rBuffer[r].begin() + i + 2, rBuffer[r].begin() + i + 2 + Element::encode((Element::CODE)rBuffer[r][i + 1]).nodes);
					edist.push_back(edist.back() + Element::encode((Element::CODE)rBuffer[r][i + 1]).nodes);
				}
			}
		}
	}

	{ // build elements
		swap(merged, clustered);
		merged.offsets.swap(clustered.offsets);
		size_t offset = 0, enodes = 0;
		for (size_t e = 0; e < merged.offsets.size(); ++e) {
			bool insert = true;
			if (Mesh::element(merged.etype[e]).dimension == dimension) {
				insert = duplicateElement.count(e) == 0;
			} else if (duplicateBoundary.count(e)) { // without duplication the parent is always local
				insert = false;
				auto p = std::lower_bound(parents.begin(), parents.end(), __parent__{ merged.offsets[e], 0, 0 });
				if (p != parents.end() && p->offset == merged.offsets[e] && p->rank == info::mpi::rank) {
					insert = true;
				}
			}
			if (insert) {
				merged.etype[offset] = merged.etype[e];
				merged.offsets[offset] = merged.offsets[e];
				for (esint n = edist[e]; n < edist[e + 1]; ++n) {
					merged.enodes[enodes++] = merged.enodes[n];
				}
				++offset;
			}
		}
		merged.offsets.resize(offset);
		merged.etype.resize(offset);
		merged.enodes.resize(enodes);
	}

	{ // build nodes
		std::unordered_map<esint, esint> dmap;
		for (size_t n = 0; n < nodes.duplication.size(); ++n) {
			dmap[nodes.duplication[n].duplicate] = g2l[nodes.duplication[n].origin];
		}
		std::unordered_map<esint, esint> rmap;
		for (size_t n = 0; n < nodes.neighbors.size(); ++n) {
			rmap[nodes.neighbors[n]] = n;
		}
		std::vector<int> neighbors(nodes.neighbors.size());

		std::vector<esint> usedNode(nodes.offsets.size() + 1);
		for (size_t e = 0, eoffset = 0; e < merged.etype.size(); eoffset += Element::encode(merged.etype[e++]).nodes) {
			PolyElement poly(merged.etype[e], merged.enodes.data() + eoffset);
			for (int n = 0; n < Element::encode(merged.etype[e]).nodes; ++n) {
				if (poly.isNode(n)) {
					usedNode[g2l[merged.enodes[n + eoffset]]] = 1;
				}
			}
		}
		utils::sizesToOffsets(usedNode);
		for (size_t e = 0, eoffset = 0; e < merged.etype.size(); eoffset += Element::encode(merged.etype[e++]).nodes) {
			PolyElement poly(merged.etype[e], merged.enodes.data() + eoffset);
			for (int n = 0; n < Element::encode(merged.etype[e]).nodes; ++n) {
				if (poly.isNode(n)) {
					merged.enodes[n + eoffset] = usedNode[g2l[merged.enodes[n + eoffset]]];
				}
			}
		}

		// we need to inform about erased nodes
		std::vector<std::vector<esint> > sRemoved(nodes.neighbors.size()), rRemoved(nodes.neighbors.size());
		for (size_t n = 0; n < nodes.offsets.size(); ++n) {
			if (usedNode[n] == usedNode[n + 1] && usedNode[g2l[nodes.offsets[n]]] == usedNode[g2l[nodes.offsets[n]] + 1]) {
				for (esint r = nodes.ranks.distribution[n]; r < nodes.ranks.distribution[n + 1]; ++r) {
					if (nodes.ranks.data[r] != info::mpi::rank) {
						sRemoved[rmap[nodes.ranks.data[r]]].push_back(nodes.offsets[n]);
					}
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sRemoved, rRemoved, nodes.neighbors)) {
			eslog::internalFailure("exchange removed nodes.\n");
		}

		std::vector<std::pair<esint, esint> > removed;
		for (size_t r = 0; r < rRemoved.size(); ++r) {
			for (size_t i = 0; i < rRemoved[r].size(); ++i) {
				auto dup = dmap.find(rRemoved[r][i]);
				if (dup != dmap.end()) {
					removed.push_back(std::pair<esint, esint>(dup->second, nodes.neighbors[r]));
				} else {
					removed.push_back(std::pair<esint, esint>(rRemoved[r][i], nodes.neighbors[r]));
				}
			}
		}
		std::sort(removed.begin(), removed.end());

		size_t offset = 0, ranks = 0, rindex = 0;
		for (size_t n = 0; n < nodes.offsets.size(); ++n) {
			if (usedNode[n] < usedNode[n + 1]) {
				while (rindex < removed.size() && removed[rindex].first < nodes.offsets[n]) { ++rindex; }
				nodes.offsets[offset] = nodes.offsets[n];
				nodes.coordinates[offset] = nodes.coordinates[n];
				size_t rdist = ranks;
				for (esint r = nodes.ranks.distribution[n]; r < nodes.ranks.distribution[n + 1]; ++r) {
					if (rindex == removed.size() || removed[rindex].first != nodes.offsets[n] || removed[rindex].second != nodes.ranks.data[r]) {
						nodes.ranks.data[ranks++] = nodes.ranks.data[r];
						if (nodes.ranks.data[r] != info::mpi::rank) {
							++neighbors[rmap[nodes.ranks.data[r]]];
						}
					}
					if (rindex < removed.size() && removed[rindex].first == nodes.offsets[n] && removed[rindex].second == nodes.ranks.data[r]) {
						++rindex;
					}
				}
				nodes.ranks.distribution[offset++] = rdist;
			}
		}
		nodes.offsets.resize(offset);
		nodes.coordinates.resize(offset);
		nodes.ranks.distribution.resize(offset + 1);
		nodes.ranks.distribution.back() = ranks;
		nodes.ranks.data.resize(ranks);
		std::vector<int> neighs; neighs.swap(nodes.neighbors);
		for (size_t n = 0; n < neighbors.size(); ++n) {
			if (neighbors[n]) {
				nodes.neighbors.push_back(neighs[n]);
			}
		}
	}

	size_t duplicates = merged.duplication.size();
	Communication::allReduce(&duplicates, nullptr, 1, MPITools::getType(duplicates).mpitype, MPI_SUM);
	eslog::info(" == DUPLICATED ELEMENTS %67lu == \n", duplicates);
}

struct __element_builder__ {
	const ivector<esint> &nodes, &dist, &owner, &neighbor;

	__element_builder__(Faces &faces, ivector<esint> &fdist, ivector<esint> &opermutation, ivector<esint> &npermutation, esint e)
	: nodes(faces.fnodes), dist(fdist), owner(faces.owner), neighbor(faces.neighbor),
	  it(owner, neighbor, opermutation, npermutation, e) { }

	struct __element__ {
		int triangles = 0, squares = 0, others = 0, enodes = 0;

		void add(const esint &nodes)
		{
			switch (nodes) {
			case 3: ++triangles; enodes += nodes + 1; break;
			case 4: ++squares; enodes += nodes + 1; break;
			default: ++others; enodes += nodes;
			}
		}

		Element::CODE code()
		{
			if (triangles == 4 && squares == 0 && others == 0) { return Element::CODE::TETRA4; }
			if (triangles == 4 && squares == 1 && others == 0) { return Element::CODE::PYRAMID5; }
			if (triangles == 2 && squares == 3 && others == 0) { return Element::CODE::PRISMA6; }
			if (triangles == 0 && squares == 6 && others == 0) { return Element::CODE::HEXA8; }
			return Element::decode(Element::CODE::POLYHEDRON, enodes + 1);
		}
	} element;

	struct __eiterator__ {
		ivector<esint>::iterator obegin, oend, nbegin, nend, olast, nlast;

		__eiterator__(const ivector<esint> &owner, const ivector<esint> &neighbor, ivector<esint> &opermutation, ivector<esint> &npermutation, esint e)
		{
			obegin = oend = std::lower_bound(opermutation.begin(), opermutation.end(), e, [&owner] (const esint &p, const esint &e) { return owner[p] < e; });
			nbegin = nend = std::lower_bound(npermutation.begin(), npermutation.end(), e, [&neighbor] (const esint &p, const esint &e) { return neighbor[p] < e; });
			olast = opermutation.end();
			nlast = npermutation.end();
		}

		void next()
		{
			obegin = oend;
			nbegin = nend;
		}

		void remove(ivector<esint>::iterator &begin, ivector<esint>::iterator &it)
		{
			std::swap(*begin , *it);
			++begin;
		}
	} it;

	void count(const esint &e)
	{
		element.triangles = element.squares = element.others = element.enodes = 0;
		while (it.oend != it.olast && owner[*it.oend] == e) {
			element.add(dist[*it.oend + 1] - dist[*it.oend]); ++it.oend;
		}
		while (it.nend != it.nlast && neighbor[*it.nend] == e) {
			element.add(dist[*it.nend + 1] - dist[*it.nend]); ++it.nend;
		}
	}

	template <int N>
	bool check(const esint *face, const std::array<esint*, N> &nodes)
	{
		for (int i = 0; i < N; ++i) {
			if (face[i] == *nodes[0] && face[(i + 1) % N] == *nodes[1]) {
				for (int j = 2; j < N; ++j) { *nodes[j] = face[(i + j) % N]; } return true;
			}
			if (face[i] == *nodes[0] && face[(i + N - 1) % N] == *nodes[1]) {
				for (int j = 2; j < N; ++j) { *nodes[j] = face[(i - j + N) % N]; } return true;
			}
		}
		return false;
	}

	template <int N>
	void finish(const std::array<esint*, N> &nodes)
	{
		for (auto ii = it.obegin; ii != it.oend; ++ii) {
			if (dist[*ii + 1] - dist[*ii] == N && check<N>(this->nodes.data() + dist[*ii], nodes)) {
				it.remove(it.obegin, ii); return;
			}
		}
		for (auto ii = it.nbegin; ii != it.nend; ++ii) {
			if (dist[*ii + 1] - dist[*ii] == N && check<N>(this->nodes.data() + dist[*ii], nodes)) {
				it.remove(it.nbegin, ii); return;
			}
		}
	}

	template <int N>
	void pushFirst(const std::array<esint*, N> &nodes)
	{
		for (auto ii = it.obegin; ii != it.oend; ++ii) {
			if (dist[*ii + 1] - dist[*ii] == N) {
				for (int i = 0; i < N; ++i) { *nodes[i] = this->nodes[dist[*ii] + i]; }
				it.remove(it.obegin, ii); return;
			}
		}
		for (auto ii = it.nbegin; ii != it.nend; ++ii) {
			if (dist[*ii + 1] - dist[*ii] == N) {
				for (int i = 0; i < N; ++i) { *nodes[i] = this->nodes[dist[*ii] + N - i - 1]; }
				it.remove(it.nbegin, ii); return;
			}
		}
	}

	void pushPolyhedron(esint* nodes)
	{
		int i = 0;
		nodes[i++] = (it.oend - it.obegin) + (it.nend - it.nbegin);
		while (it.obegin != it.oend) {
			if (dist[*it.obegin + 1] - dist[*it.obegin] == 3 || dist[*it.obegin + 1] - dist[*it.obegin] == 4) {
				nodes[i++] = dist[*it.obegin + 1] - dist[*it.obegin];
			}
			for (auto n = dist[*it.obegin]; n < dist[*it.obegin + 1]; ++n) {
				nodes[i++] = this->nodes[n];
			}
			++it.obegin;
		}
		while (it.nbegin != it.nend) {
			int header = 0;
			if (dist[*it.nbegin + 1] - dist[*it.nbegin] == 3 || dist[*it.nbegin + 1] - dist[*it.nbegin] == 4) {
				nodes[i++] = dist[*it.nbegin + 1] - dist[*it.nbegin];
			} else {
				header = 1;
				nodes[i++] = this->nodes[dist[*it.nbegin]];
			}

			for (auto n = dist[*it.nbegin + 1]; n > dist[*it.nbegin] + header; --n) {
				nodes[i++] = this->nodes[n - 1];
			}
			++it.nbegin;
		}
	}
};

void buildElementsFromFaces(Faces &faces, Elements &elements)
{
	esint poly[2] = { 0, 0 };
	if (faces.ftype.size()) {
		ivector<esint> opermutation(faces.owner.size()), npermutation(faces.neighbor.size());
		std::iota(opermutation.begin(), opermutation.end(), 0);
		std::sort(opermutation.begin(), opermutation.end(), [&] (const esint &i, const esint &j) { return faces.owner[i] != faces.owner[j] ? faces.owner[i] < faces.owner[j] : i < j; });
		std::iota(npermutation.begin(), npermutation.end(), 0);
		std::sort(npermutation.begin(), npermutation.end(), [&] (const esint &i, const esint &j) { return faces.neighbor[i] != faces.neighbor[j] ? faces.neighbor[i] < faces.neighbor[j] : i < j; });

		esint min = faces.owner[opermutation.front()], max = std::max(faces.owner[opermutation.back()], faces.neighbor[npermutation.back()]);
		for (auto n = npermutation.begin(); n != npermutation.end(); ++n) {
			if (faces.neighbor[*n] != -1) {
				min = std::min(min, faces.neighbor[*n]);
			}
		}

		ivector<esint> fdist(faces.ftype.size() + 1);
		fdist[0] = 0;
		for (size_t e = 0; e < faces.ftype.size(); ++e) {
			if (Element::encode(faces.ftype[e]).code == Element::CODE::POLYGON) {
				++poly[0];
			}
			fdist[e + 1] = fdist[e] + Element::encode(faces.ftype[e]).nodes;
		}

		elements.etype.reserve(max - min + 1);
		elements.enodes.reserve(faces.fnodes.size());
		__element_builder__ builder(faces, fdist, opermutation, npermutation, min);
		for (esint e = min; e <= max; ++e, builder.it.next()) {
			builder.count(e);
			elements.etype.push_back(builder.element.code());
			elements.enodes.resize(elements.enodes.size() + Element::encode(elements.etype.back()).nodes, -1);
			auto enodes = elements.enodes.data() + elements.enodes.size() - Element::encode(elements.etype.back()).nodes;
			switch (elements.etype.back()) {
			case Element::CODE::TETRA4:
				builder.pushFirst<3>({ enodes, enodes + 1, enodes + 2 });
				builder.finish<3>({ enodes, enodes + 1, enodes + 3 });
				break;
			case Element::CODE::PYRAMID5:
				builder.pushFirst<4>({ enodes, enodes + 1, enodes + 2, enodes + 3 });
				builder.finish<3>({ enodes, enodes + 1, enodes + 4 });
				break;
			case Element::CODE::PRISMA6:
				builder.pushFirst<3>({ enodes + 3, enodes + 4, enodes + 5 });
				builder.finish<4>({ enodes + 4, enodes + 3, enodes + 0, enodes + 1 });
				builder.finish<3>({ enodes + 0, enodes + 1, enodes + 2 });
				break;
			case Element::CODE::HEXA8:
				builder.pushFirst<4>({ enodes, enodes + 1, enodes + 5, enodes + 4 });
				builder.finish<4>({ enodes + 5, enodes + 1, enodes + 2, enodes + 6 });
				builder.finish<4>({ enodes + 0, enodes + 4, enodes + 7, enodes + 3 });
				break;
			default:
				++poly[1];
				builder.pushPolyhedron(enodes);
			}
		}
	}

	Communication::allReduce(poly, nullptr, 2, MPITools::getType<esint>().mpitype, MPI_SUM);
	eslog::info(" == TOTAL NUMBER OF POLYGONS %62d == \n", poly[0]);
	eslog::info(" == TOTAL NUMBER OF POLYHEDRONS %59d == \n", poly[1]);
}

void rotateNormalsOut(LinkedNodes &nodes, MergedElements &elements)
{
	std::vector<esint> edistribution = tarray<esint>::distribute(info::env::OMP_NUM_THREADS, elements.etype.size());

	ivector<esint> edist(elements.etype.size() + 1);
	edist[0] = 0;
	for (size_t e = 0; e < elements.etype.size(); ++e) {
		edist[e + 1] = edist[e] + Element::encode(elements.etype[e]).nodes;
	}

	auto dist = [] (const Point &p, const Point &v1, const Point &v2, const Point &v3) {
		return (p - v1) * Point::cross((v2 - v1).normalize(), (v3 - v1).normalize());
	};

	auto distance = [&nodes,&dist] (esint t1, esint t2, esint t3, esint p) { return dist(nodes.coordinates[p], nodes.coordinates[t1], nodes.coordinates[t2], nodes.coordinates[t3]); };
	auto distancep = [&nodes,&dist] (esint t1, esint t2, esint t3, const Point &p) { return dist(p, nodes.coordinates[t1], nodes.coordinates[t2], nodes.coordinates[t3]); };

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		for (esint e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			esint *enodes = elements.enodes.data() + edist[e];
			switch (Element::encode(elements.etype[e]).code) {
			case Element::CODE::TETRA4:
				if (distance(enodes[0], enodes[1], enodes[2], enodes[3]) < 0) {
					std::swap(enodes[1], enodes[2]);
				}
				break;
			case Element::CODE::PYRAMID5:
				if (distance(enodes[0], enodes[1], enodes[2], enodes[4]) < 0) {
					std::swap(enodes[0], enodes[2]);
				}
				break;
			case Element::CODE::PRISMA6:
				if (distance(enodes[0], enodes[1], enodes[2], enodes[3]) < 0) {
					std::swap(enodes[0], enodes[3]);
					std::swap(enodes[1], enodes[4]);
					std::swap(enodes[2], enodes[5]);
				}
				break;
			case Element::CODE::HEXA8:
				if (distance(enodes[0], enodes[1], enodes[2], enodes[4]) < 0) {
					std::swap(enodes[0], enodes[4]);
					std::swap(enodes[1], enodes[5]);
					std::swap(enodes[2], enodes[6]);
					std::swap(enodes[3], enodes[7]);
				}
				break;
			case Element::CODE::POLYHEDRON: {
				Point center;
				PolyElement poly(elements.etype[e], enodes);
				for (int n = 0; n < poly.size; ++n) {
					if(poly.isNode(n)) {
						center += nodes.coordinates[enodes[n]];
					}
				}
				center /= poly.nodes;
				for (int p = 0, n = 1; p < enodes[0]; ++p, n += enodes[n] + 1) {
					if (enodes[n] == 3) {
						if (distancep(enodes[n + 1], enodes[n + 2], enodes[n + 3], center) < 0) {
							std::swap(enodes[n + 2], enodes[n + 3]);
						}
					} else {
						Point pcenter;
						for (esint pn = 0; pn < enodes[n]; ++pn) {
							pcenter += nodes.coordinates[enodes[n + pn + 1]];
						}
						pcenter /= enodes[n];
						if (dist(pcenter, nodes.coordinates[enodes[n + 1]], nodes.coordinates[enodes[n + 2]], center) < 0) {
							for (int nn = 0; nn < enodes[n] / 2; ++nn) {
								std::swap(enodes[n + 1 + nn], enodes[n + enodes[n] - 1 - nn]);
							}
						}
					}
				}
			} break;
			default:
				eslog::error("not implemented normal rotation for the loaded element type: %d\n", (int)elements.etype[e]);
			}
		}
	}
}

}
}
