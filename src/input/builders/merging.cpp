
#include "builder.utils.h"

#include "basis/containers/tarray.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/structures/kdtree.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/envinfo.h"
#include "wrappers/mpi/communication.h"

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

static void localSearch(MergedNodes *merged)
{
	treeSearch(merged->coordinates, [&] (esint origin, esint target) {
		merged->duplication.push_back({ target, origin });
	});
	std::sort(merged->duplication.begin(), merged->duplication.end());

	// set origin to the last occurrence -> it simplify linking step
	for (auto begin = merged->duplication.begin(), end = begin; begin != merged->duplication.end(); begin = end) {
		while (end != merged->duplication.end() && begin->origin == end->origin) { ++end; }
		esint max = std::max((end - 1)->origin, (end - 1)->duplicate);
		for (auto it = begin; it != end; ++it) {
			if (it->duplicate == max) {
				std::swap(it->duplicate, it->origin);
			} else {
				it->origin = max;
			}
		}
	}
	std::sort(merged->duplication.begin(), merged->duplication.end());
}

void searchDuplicatedNodes(TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, TemporalSequentialMesh<MergedNodes, ClusteredElements> &merged)
{
	if (info::mpi::size != 1) {
		eslog::internalFailure("usage of a sequential method during a parallel run.\n");
	}
	merged.nodes->coordinates.swap(clustered.nodes->coordinates);
	merged.nodes->offsets.swap(clustered.nodes->offsets);
	localSearch(merged.nodes);

	delete merged.elements;
	merged.elements = clustered.elements;
	clustered.elements = nullptr;
	clustered.clear();
}

void searchDuplicatedNodes(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, const std::vector<int> &sfcNeighbors, TemporalMesh<ClusteredNodes, ClusteredElements> &clustered, TemporalMesh<MergedNodes, ClusteredElements> &merged)
{
	std::vector<esint> cdistribution = tarray<esint>::distribute(info::env::OMP_NUM_THREADS, clustered.nodes->offsets.size());

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
				size_t origin = sfc.getBucket(clustered.nodes->coordinates[n][d], d);
				if (sfc.getBucket(clustered.nodes->coordinates[n][d] - eps, d) < origin) {
					hit[2 * d + 0] = -1;
				}
				if (sfc.getBucket(clustered.nodes->coordinates[n][d] + eps, d) > origin) {
					hit[2 * d + 1] = +1;
				}
			}

			neighs.clear();
			for (int x = hit[0]; x <= hit[1]; ++x) {
				for (int y = hit[2]; y <= hit[3]; ++y) {
					for (int z = hit[4]; z <= hit[5]; ++z) {
						if (x || y || z) {
							size_t b = sfc.getBucket(clustered.nodes->coordinates[n] + Point(x * eps, y * eps, z * eps));
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
	std::vector<_Point<esfloat>, initless_allocator<_Point<esfloat> > > coordinates;
	std::vector<std::vector<esint> > sBuffer(sfcNeighbors.size()), rBuffer(sfcNeighbors.size());
	{ // build send buffer
		auto ni = sfcNeighbors.begin();
		for (auto begin = toneighs[0].begin(), end = begin; end != toneighs[0].end(); begin = end++) {
			while (*ni < begin->target) { ++ni; }
			while (end != toneighs[0].end() && begin->target == end->target) { ++end; }

			for (auto n = begin; n != end; ++n) {
				offsets.push_back(n->offset);
				coordinates.push_back(clustered.nodes->coordinates[n->offset]);
			}
			esint nodes = end - begin;
			sBuffer[*ni].resize(1 + nodes + utils::reinterpret_size<esint, _Point<esfloat> >(nodes));
			sBuffer[*ni][0] = nodes;
			memcpy(sBuffer[*ni].data() + 1, offsets.data() + offsets.size() - nodes, nodes * sizeof(esint));
			memcpy(sBuffer[*ni].data() + 1 + nodes, reinterpret_cast<esint*>(coordinates.data() + coordinates.size() - nodes), nodes * sizeof(_Point<esfloat>));
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

	size_t size = clustered.nodes->offsets.size();
	ClusteredNodes *from, *to;
	if (toerase.size() >= toinsert.size()) { // reuse memory
		clustered.nodes->offsets.swap(merged.nodes->offsets);
		clustered.nodes->coordinates.swap(merged.nodes->coordinates);
		from = to = merged.nodes;
	} else {
		merged.nodes->offsets.resize(size + toinsert.size() - toerase.size());
		merged.nodes->coordinates.resize(size + toinsert.size() - toerase.size());
		from = clustered.nodes;
		to = merged.nodes;
	}

	if (toerase.size()) {
		for (size_t n = toerase.front(), last = toerase.front(), e = 0; n < from->offsets.size(); ++n) {
			if (toerase.size() <= e || n != toerase[e]) {
				to->offsets[last] = from->offsets[n];
				to->coordinates[last] = from->coordinates[n];
				++last;
			} else {
				++e;
			}
		}
		size -= toerase.size();
		to->coordinates.resize(size);
		to->offsets.resize(size);
	}

	if (toinsert.size()) {
		to->coordinates.resize(size + toinsert.size());
		to->offsets.resize(size + toinsert.size());
		for (size_t i = size - toinsert.size(), j = 0; j < toinsert.size(); ++i, ++j) {
			to->offsets[i] = offsets[toinsert[j]];
			to->coordinates[i] = coordinates[toinsert[j]];
		}
	}

	localSearch(merged.nodes);

	delete merged.elements;
	merged.elements = clustered.elements;
	clustered.elements = nullptr;
	clustered.clear();
}

void searchDuplicatedElements(TemporalSequentialMesh<MergedNodes, ClusteredElements> &merged, TemporalSequentialMesh<MergedNodes, MergedElements> &prepared, int meshDimension)
{
	{ // build nodes
		ivector<esint> usedNodes(merged.nodes->offsets.size());
		std::iota(usedNodes.begin(), usedNodes.end(), 0);
		for (size_t n = 0; n < merged.nodes->duplication.size(); ++n) {
			usedNodes[merged.nodes->duplication[n].duplicate] = merged.nodes->duplication[n].origin;
		}

		prepared.nodes->offsets.swap(merged.nodes->offsets);
		prepared.nodes->coordinates.swap(merged.nodes->coordinates);

		for (size_t n = 0, offset = 0; n < usedNodes.size(); ++n) {
			if (usedNodes[n] == (esint)n) {
				usedNodes[n] = offset++;
			}
		}

		size_t offset = 0;
		for (size_t n = 0; n < usedNodes.size(); ++n) {
			if (usedNodes[n] <= (esint)n) {
				prepared.nodes->offsets[offset] = prepared.nodes->offsets[n];
				prepared.nodes->coordinates[offset] = prepared.nodes->coordinates[n];
				++offset;
			} else {
				usedNodes[n] = usedNodes[usedNodes[n]];
			}
		}
		prepared.nodes->offsets.resize(offset);
		prepared.nodes->coordinates.resize(offset);

		for (size_t n = 0; n < merged.elements->enodes.size(); ++n) {
			merged.elements->enodes[n] = usedNodes[merged.elements->enodes[n]];
		}
	}

	std::unordered_set<esint> duplicateElement;

	auto getMinimal = [&] (esint nsize, esint *nodes) {
		esint min = nodes[0];
		for (int n = 1; n < nsize; ++n) {
			esint current = nodes[n];
			if (prepared.nodes->coordinates[current].x < prepared.nodes->coordinates[min].x) {
				min = current;
			} else if (prepared.nodes->coordinates[current].x == prepared.nodes->coordinates[min].x && current < min) {
				min = current;
			}
		}
		return min;
	};

	std::vector<esint> mapDist(prepared.nodes->offsets.size() + 1);
	ivector<esint> mapData;

	for (size_t e = 0; e < merged.elements->offsets.size(); ++e) {
		if (Mesh::element(merged.elements->etype[e]).dimension == meshDimension) {
			++mapDist[getMinimal(merged.elements->edist[e + 1] - merged.elements->edist[e], merged.elements->enodes.data() + merged.elements->edist[e])];
		}
	}
	utils::sizesToOffsets(mapDist);
	mapData.resize(mapDist.back());

	std::vector<esint> _mapDist = mapDist;
	for (size_t e = 0; e < merged.elements->offsets.size(); ++e) {
		if (Mesh::element(merged.elements->etype[e]).dimension == meshDimension) {
			mapData[_mapDist[getMinimal(merged.elements->edist[e + 1] - merged.elements->edist[e], merged.elements->enodes.data() + merged.elements->edist[e])]++] = e;
		}
	}
	utils::clearVector(_mapDist);

	std::vector<esint> _checkBuffer;
	_checkBuffer.reserve(2 * 20);
	auto areSame = [&] (esint size1, esint *nodes1, esint size2, esint *nodes2) {
		if (size1 == size2) {
			_checkBuffer.clear();
			for (esint n = 0; n < size1; ++n) { _checkBuffer.push_back(nodes1[n]); }
			for (esint n = 0; n < size2; ++n) { _checkBuffer.push_back(nodes2[n]); }
			std::sort(_checkBuffer.begin(), _checkBuffer.begin() + size1);
			std::sort(_checkBuffer.begin() + size1, _checkBuffer.end());
			return memcmp(_checkBuffer.data(), _checkBuffer.data() + size1, sizeof(esint) * size1) == 0;
		}
		return false;
	};
	for (size_t e1 = 0; e1 < merged.elements->offsets.size(); ++e1) {
		if (Mesh::element(merged.elements->etype[e1]).dimension == meshDimension && duplicateElement.count(e1) == 0) {
			esint min = getMinimal(merged.elements->edist[e1 + 1] - merged.elements->edist[e1], merged.elements->enodes.data() + merged.elements->edist[e1]);
			for (esint n = mapDist[min]; n < mapDist[min + 1]; ++n) {
				size_t e2 = mapData[n];
				if (e1 < e2) { // if e2 is lower, the pair was already processed
					if (Mesh::element(merged.elements->etype[e2]).dimension == meshDimension && areSame(merged.elements->edist[e1 + 1] - merged.elements->edist[e1], merged.elements->enodes.data() + merged.elements->edist[e1], merged.elements->edist[e2 + 1] - merged.elements->edist[e2], merged.elements->enodes.data() + merged.elements->edist[e2])) {
						prepared.elements->duplication.push_back({ merged.elements->offsets[e1], merged.elements->offsets[e2] });
						duplicateElement.insert(e2);
					}
				}
			}
		}
	}
	std::sort(prepared.elements->duplication.begin(), prepared.elements->duplication.end());

	{ // build elements
		prepared.elements->offsets.swap(merged.elements->offsets);
		prepared.elements->etype.swap(merged.elements->etype);
		prepared.elements->enodes.swap(merged.elements->enodes);
		prepared.elements->edist.swap(merged.elements->edist);
		size_t offset = 0, enodes = 0;
		for (size_t e = 0; e < prepared.elements->offsets.size(); ++e) {
			if (duplicateElement.count(e) == 0) {
				prepared.elements->etype[offset] = prepared.elements->etype[e];
				prepared.elements->offsets[offset] = prepared.elements->offsets[e];
				for (esint n = prepared.elements->edist[e]; n < prepared.elements->edist[e + 1]; ++n) {
					prepared.elements->enodes[enodes++] = prepared.elements->enodes[n];
				}
				prepared.elements->edist[offset] = enodes - (prepared.elements->edist[e + 1] - prepared.elements->edist[e]);
				++offset;
			}
		}
		prepared.elements->offsets.resize(offset);
		prepared.elements->etype.resize(offset);
		prepared.elements->edist.resize(offset + 1);
		prepared.elements->edist.back() = enodes;
		prepared.elements->enodes.resize(enodes);
	}

	merged.clear();
}

void searchParentAndDuplicatedElements(TemporalMesh<LinkedNodes, ClusteredElements> &linked, TemporalMesh<LinkedNodes, MergedElements> &prepared, int meshDimension)
{
	// duplicate boundary is adept to duplication only
	std::unordered_set<esint> duplicateBoundary, duplicateElement;

	// eoffset, etype, enodes; .... + project all enodes
	std::vector<std::vector<esint> > sBuffer(linked.nodes->neighbors.size()), rBuffer(linked.nodes->neighbors.size());
	if (linked.nodes->neighbors.size()) { // duplicated elements have to have all nodes held by other rank
		std::vector<int> neighbors;
		for (size_t e = 0; e < linked.elements->offsets.size(); ++e) {
			neighbors.assign(linked.nodes->neighbors.begin(), linked.nodes->neighbors.end());
			for (esint n = linked.elements->edist[e]; n < linked.elements->edist[e + 1]; ++n) {
				esint local = linked.nodes->g2l[linked.elements->enodes[n]];
				size_t intersection = 0, current = 0;
				for (int r = linked.nodes->rankDistribution[local]; r < linked.nodes->rankDistribution[local + 1]; ++r) {
					while (current < neighbors.size() && neighbors[current] < linked.nodes->rankData[r]) { ++current; }
					if (current < neighbors.size() && neighbors[current] == linked.nodes->rankData[r]) {
						neighbors[intersection++] = neighbors[current++];
					}
				}
				neighbors.resize(intersection);
			}
			esint roffset = 0;
			for (auto r = neighbors.begin(); r != neighbors.end(); ++r) {
				while (linked.nodes->neighbors[roffset] < *r) { ++roffset; }
				sBuffer[roffset].push_back(linked.elements->offsets[e]);
				sBuffer[roffset].push_back((esint)linked.elements->etype[e]);
				for (esint n = linked.elements->edist[e]; n < linked.elements->edist[e + 1]; ++n) {
					sBuffer[roffset].push_back(linked.elements->enodes[n]);
				}
			}
			if (neighbors.size() && Mesh::element(linked.elements->etype[e]).dimension < meshDimension) {
				duplicateBoundary.insert(e);
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, linked.nodes->neighbors)) {
		eslog::internalFailure("cannot exchange duplicated elements.\n");
	}
	utils::clearVector(sBuffer);

	auto getMinimal = [&] (esint nsize, esint *nodes) {
		esint min = linked.nodes->g2l[nodes[0]];
		for (int n = 1; n < nsize; ++n) {
			esint current = linked.nodes->g2l[nodes[n]];
			if (linked.nodes->coordinates[current].x < linked.nodes->coordinates[min].x) {
				min = current;
			} else if (linked.nodes->coordinates[current].x == linked.nodes->coordinates[min].x && current < min) {
				min = current;
			}
		}
		return min;
	};

	std::vector<esint> mapDist(linked.nodes->offsets.size() + 1);
	ivector<esint> mapData;

	for (size_t e = 0; e < linked.elements->offsets.size(); ++e) {
		if (Mesh::element(linked.elements->etype[e]).dimension == meshDimension) {
			++mapDist[getMinimal(linked.elements->edist[e + 1] - linked.elements->edist[e], linked.elements->enodes.data() + linked.elements->edist[e])];
		} else {
			if (duplicateBoundary.count(e)) { // we have to search a parent for each sent elements
				for (esint n = linked.elements->edist[e]; n < linked.elements->edist[e + 1]; ++n) {
					++mapDist[linked.nodes->g2l[linked.elements->enodes[n]]];
				}
			}
		}
	}
	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Mesh::element(rBuffer[r][i + 1]).nodes) {
			if (Mesh::element(rBuffer[r][i + 1]).dimension == meshDimension) {
				++mapDist[getMinimal(Mesh::element(rBuffer[r][i + 1]).nodes, rBuffer[r].data() + i + 2)];
			} else {
				for (int n = 0; n < Mesh::element(rBuffer[r][i + 1]).nodes; ++n) {
					++mapDist[linked.nodes->g2l[rBuffer[r][i + 2 + n]]];
				}
			}
		}
	}
	utils::sizesToOffsets(mapDist);
	mapData.resize(mapDist.back());

	std::vector<std::pair<esint, esint> > recvMap;
	std::vector<esint> _mapDist = mapDist;

	for (size_t e = 0; e < linked.elements->offsets.size(); ++e) {
		if (Mesh::element(linked.elements->etype[e]).dimension == meshDimension) {
			mapData[_mapDist[getMinimal(linked.elements->edist[e + 1] - linked.elements->edist[e], linked.elements->enodes.data() + linked.elements->edist[e])]++] = e;
		} else {
			if (duplicateBoundary.count(e)) {
				for (esint n = linked.elements->edist[e]; n < linked.elements->edist[e + 1]; ++n) {
					mapData[_mapDist[linked.nodes->g2l[linked.elements->enodes[n]]]++] = e;
				}
			}
		}
	}
	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Mesh::element(rBuffer[r][i + 1]).nodes) {
			if (Mesh::element(rBuffer[r][i + 1]).dimension == meshDimension) {
				mapData[_mapDist[getMinimal(Mesh::element(rBuffer[r][i + 1]).nodes, rBuffer[r].data() + i + 2)]++] = linked.elements->offsets.size() + recvMap.size();
			} else {
				for (int n = 0; n < Mesh::element(rBuffer[r][i + 1]).nodes; ++n) {
					mapData[_mapDist[linked.nodes->g2l[rBuffer[r][i + 2 + n]]]++] = linked.elements->offsets.size() + recvMap.size();
				}
			}
			recvMap.push_back({ (esint)r, (esint)i });
		}
	}
	utils::clearVector(_mapDist);

	ivector<esint> shared, tmpShared;
	auto merge = [&] (esint nsize, esint *nodes) {
		esint min = linked.nodes->g2l[nodes[0]];
		shared.assign(mapData.begin() + mapDist[min], mapData.begin() + mapDist[min + 1]);
		for (esint n = 1; n < nsize; ++n) {
			esint current = linked.nodes->g2l[nodes[n]];
			tmpShared.swap(shared);
			shared.resize(tmpShared.size() + mapDist[current + 1] - mapDist[current]);
			std::merge(mapData.begin() + mapDist[current], mapData.begin() + mapDist[current + 1], tmpShared.begin(), tmpShared.end(), shared.begin());
			if (linked.nodes->coordinates[current].x < linked.nodes->coordinates[min].x) {
				min = current;
			} else if (linked.nodes->coordinates[current].x == linked.nodes->coordinates[min].x && current < min) {
				min = current;
			}
		}
		return min;
	};

	struct __parent__ {
		esint offset, parent, rank;
		bool operator<(const __parent__ &other) const { return offset == other.offset ? parent < other.parent : offset < other.offset; }
	};

	std::vector<__parent__> parents;
	std::vector<std::vector<__parent__> > rParents(linked.nodes->neighbors.size());
	std::vector<esint> _checkBuffer;
	_checkBuffer.reserve(2 * 20);
	auto areSame = [&] (esint size1, esint *nodes1, esint size2, esint *nodes2) {
		if (size1 == size2) {
			_checkBuffer.clear();
			for (esint n = 0; n < size1; ++n) { _checkBuffer.push_back(linked.nodes->g2l[nodes1[n]]); }
			for (esint n = 0; n < size2; ++n) { _checkBuffer.push_back(linked.nodes->g2l[nodes2[n]]); }
			std::sort(_checkBuffer.begin(), _checkBuffer.begin() + size1);
			std::sort(_checkBuffer.begin() + size1, _checkBuffer.end());
			return memcmp(_checkBuffer.data(), _checkBuffer.data() + size1, sizeof(esint) * size1) == 0;
		}
		return false;
	};
	for (size_t e1 = 0; e1 < linked.elements->offsets.size(); ++e1) {
		if (Mesh::element(linked.elements->etype[e1]).dimension == meshDimension && duplicateElement.count(e1) == 0) {
			esint min = merge(linked.elements->edist[e1 + 1] - linked.elements->edist[e1], linked.elements->enodes.data() + linked.elements->edist[e1]);
			for (esint n = mapDist[min]; n < mapDist[min + 1]; ++n) {
				size_t e2 = mapData[n];
				if (e1 < e2) { // if e2 is lower, the pair was already processed
					if (e2 < linked.elements->offsets.size()) { // local
						if (Mesh::element(linked.elements->etype[e2]).dimension == meshDimension && areSame(linked.elements->edist[e1 + 1] - linked.elements->edist[e1], linked.elements->enodes.data() + linked.elements->edist[e1], linked.elements->edist[e2 + 1] - linked.elements->edist[e2], linked.elements->enodes.data() + linked.elements->edist[e2])) {
							prepared.elements->duplication.push_back({ linked.elements->offsets[e1], linked.elements->offsets[e2] });
							duplicateElement.insert(e2);
						}
					} else { // from neighbors
						const auto &r = recvMap[e2 - linked.elements->offsets.size()];
						if (Mesh::element(rBuffer[r.first][r.second + 1]).dimension == meshDimension) {
							if (areSame(linked.elements->edist[e1 + 1] - linked.elements->edist[e1], linked.elements->enodes.data() + linked.elements->edist[e1], Mesh::element(rBuffer[r.first][r.second + 1]).nodes, rBuffer[r.first].data() + r.second + 2)) {
								if (linked.elements->offsets[e1] < rBuffer[r.first][r.second]) {
									prepared.elements->duplication.push_back({ linked.elements->offsets[e1], rBuffer[r.first][r.second] });
									// element is removed on the neighboring process
								} else {
									prepared.elements->duplication.push_back({ rBuffer[r.first][r.second], linked.elements->offsets[e1] });
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
					if (e2 < linked.elements->offsets.size()) {
						if (Mesh::element(linked.elements->etype[e2]).nodes == end - begin) {
							parents.push_back({ linked.elements->offsets[e2], linked.elements->offsets[e1], info::mpi::rank });
						}
					} else {
						auto &roffset = recvMap[e2 - linked.elements->offsets.size()];
						if (Mesh::element(rBuffer[roffset.first][roffset.second + 1]).nodes == end - begin) {
							parents.push_back({ rBuffer[roffset.first][roffset.second], linked.elements->offsets[e1], info::mpi::rank });
						}
					}
				}
			}
		}
	}

	if (!Communication::exchangeUnknownSize(parents, rParents, linked.nodes->neighbors)) {
		eslog::internalFailure("exchange parent elements.\n");
	}
	for (size_t r = 0; r < rParents.size(); ++r) {
		parents.insert(parents.end(), rParents[r].begin(), rParents[r].end());
	}
	std::sort(parents.begin(), parents.end());
	std::sort(prepared.elements->duplication.begin(), prepared.elements->duplication.end());

	// linked.elements should have enough capacity to insert some other elements (see clusterization)
	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Mesh::element(rBuffer[r][i + 1]).nodes) {
			if (Mesh::element(rBuffer[r][i + 1]).dimension < meshDimension) { // only boundary elements can be inserted
				auto p = std::lower_bound(parents.begin(), parents.end(), __parent__{ rBuffer[r][i], 0, 0 });
				if (p != parents.end() && p->offset == rBuffer[r][i] && p->rank == info::mpi::rank) {
					linked.elements->offsets.push_back(rBuffer[r][i]);
					linked.elements->etype.push_back((Element::CODE)rBuffer[r][i + 1]);
					linked.elements->enodes.insert(linked.elements->enodes.end(), rBuffer[r].begin() + i + 2, rBuffer[r].begin() + i + 2 + Mesh::element(rBuffer[r][i + 1]).nodes);
					linked.elements->edist.push_back(linked.elements->enodes.size());
				}
			}
		}
	}

	{ // build elements
		prepared.elements->offsets.swap(linked.elements->offsets);
		prepared.elements->etype.swap(linked.elements->etype);
		prepared.elements->enodes.swap(linked.elements->enodes);
		prepared.elements->edist.swap(linked.elements->edist);
		size_t offset = 0, enodes = 0;
		for (size_t e = 0; e < prepared.elements->offsets.size(); ++e) {
			bool insert = true;
			if (Mesh::element(prepared.elements->etype[e]).dimension == meshDimension) {
				insert = duplicateElement.count(e) == 0;
			} else {
				insert = false;
				auto p = std::lower_bound(parents.begin(), parents.end(), __parent__{ prepared.elements->offsets[e], 0, 0 });
				if (p != parents.end() && p->offset == prepared.elements->offsets[e] && p->rank == info::mpi::rank) {
					insert = true;
				}
			}
			if (insert) {
				prepared.elements->etype[offset] = prepared.elements->etype[e];
				prepared.elements->offsets[offset] = prepared.elements->offsets[e];
				for (esint n = prepared.elements->edist[e]; n < prepared.elements->edist[e + 1]; ++n) {
					prepared.elements->enodes[enodes++] = prepared.elements->enodes[n];
				}
				prepared.elements->edist[offset] = enodes - (prepared.elements->edist[e + 1] - prepared.elements->edist[e]);
				++offset;
			}
		}
		prepared.elements->offsets.resize(offset);
		prepared.elements->etype.resize(offset);
		prepared.elements->edist.resize(offset + 1);
		prepared.elements->edist.back() = enodes;
		prepared.elements->enodes.resize(enodes);
	}

	{ // build nodes
		std::unordered_map<esint, esint> dmap;
		for (size_t n = 0; n < linked.nodes->duplication.size(); ++n) {
			dmap[linked.nodes->duplication[n].duplicate] = linked.nodes->g2l[linked.nodes->duplication[n].origin];
		}
		std::unordered_map<esint, esint> rmap;
		for (size_t n = 0; n < linked.nodes->neighbors.size(); ++n) {
			rmap[linked.nodes->neighbors[n]] = n;
		}
		std::vector<int> neighbors(linked.nodes->neighbors.size());

		prepared.nodes->offsets.swap(linked.nodes->offsets);
		prepared.nodes->coordinates.swap(linked.nodes->coordinates);
		prepared.nodes->rankDistribution.swap(linked.nodes->rankDistribution);
		prepared.nodes->rankData.swap(linked.nodes->rankData);

		std::vector<esint> usedNode(prepared.nodes->offsets.size() + 1);
		for (size_t n = 0; n < prepared.elements->enodes.size(); ++n) {
			usedNode[linked.nodes->g2l[prepared.elements->enodes[n]]] = 1;
		}
		utils::sizesToOffsets(usedNode);
		for (size_t n = 0; n < prepared.elements->enodes.size(); ++n) {
			prepared.elements->enodes[n] = usedNode[linked.nodes->g2l[prepared.elements->enodes[n]]];
		}

		// we need to inform about erased nodes
		std::vector<std::vector<esint> > sRemoved(linked.nodes->neighbors.size()), rRemoved(linked.nodes->neighbors.size());
		for (size_t n = 0; n < prepared.nodes->offsets.size(); ++n) {
			if (usedNode[n] == usedNode[n + 1] && usedNode[linked.nodes->g2l[prepared.nodes->offsets[n]]] == usedNode[linked.nodes->g2l[prepared.nodes->offsets[n]] + 1]) {
				for (esint r = prepared.nodes->rankDistribution[n]; r < prepared.nodes->rankDistribution[n + 1]; ++r) {
					if (prepared.nodes->rankData[r] != info::mpi::rank) {
						sRemoved[rmap[prepared.nodes->rankData[r]]].push_back(prepared.nodes->offsets[n]);
					}
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sRemoved, rRemoved, linked.nodes->neighbors)) {
			eslog::internalFailure("exchange removed nodes.\n");
		}

		std::vector<std::pair<esint, esint> > removed;
		for (size_t r = 0; r < rRemoved.size(); ++r) {
			for (size_t i = 0; i < rRemoved[r].size(); ++i) {
				auto dup = dmap.find(rRemoved[r][i]);
				if (dup != dmap.end()) {
					removed.push_back(std::pair<esint, esint>(dup->second, linked.nodes->neighbors[r]));
				} else {
					removed.push_back(std::pair<esint, esint>(rRemoved[r][i], linked.nodes->neighbors[r]));
				}
			}
		}
		std::sort(removed.begin(), removed.end());

		size_t offset = 0, ranks = 0, rindex = 0;
		for (size_t n = 0; n < prepared.nodes->offsets.size(); ++n) {
			if (usedNode[n] < usedNode[n + 1]) {
				while (rindex < removed.size() && removed[rindex].first < prepared.nodes->offsets[n]) { ++rindex; }
				prepared.nodes->offsets[offset] = prepared.nodes->offsets[n];
				prepared.nodes->coordinates[offset] = prepared.nodes->coordinates[n];
				size_t rdist = ranks;
				for (esint r = prepared.nodes->rankDistribution[n]; r < prepared.nodes->rankDistribution[n + 1]; ++r) {
					if (rindex == removed.size() || removed[rindex].first != prepared.nodes->offsets[n] || removed[rindex].second != prepared.nodes->rankData[r]) {
						prepared.nodes->rankData[ranks++] = prepared.nodes->rankData[r];
						if (prepared.nodes->rankData[r] != info::mpi::rank) {
							++neighbors[rmap[prepared.nodes->rankData[r]]];
						}
					}
					if (rindex < removed.size() && removed[rindex].first == prepared.nodes->offsets[n] && removed[rindex].second == prepared.nodes->rankData[r]) {
						++rindex;
					}
				}
				prepared.nodes->rankDistribution[offset++] = rdist;
			}
		}
		prepared.nodes->offsets.resize(offset);
		prepared.nodes->coordinates.resize(offset);
		prepared.nodes->rankDistribution.resize(offset + 1);
		prepared.nodes->rankDistribution.back() = ranks;
		prepared.nodes->rankData.resize(ranks);
		for (size_t n = 0; n < neighbors.size(); ++n) {
			if (neighbors[n]) {
				prepared.nodes->neighbors.push_back(linked.nodes->neighbors[n]);
			}
		}
	}
	linked.clear();
}

void buildElementsFromFaces()
{

}

}
}
