
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

static void searchDuplicateNodes(std::vector<_Point<esfloat>, initless_allocator<_Point<esfloat> > > &coordinates, std::function<void(esint origin, esint target)> merge)
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

void searchDuplicatedNodes(const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered)
{
	std::vector<esint> cdistribution = tarray<esint>::distribute(info::env::OMP_NUM_THREADS, clustered.noffsets.size());

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
			for (int d = 0; d < clustered.dimension; ++d) {
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
							int nn = std::lower_bound(clustered.splitters.begin(), clustered.splitters.end(), b + 1) - clustered.splitters.begin() - 1;
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
	std::vector<std::vector<esint> > sBuffer(clustered.neighbors.size()), rBuffer(clustered.neighbors.size());
	{ // build send buffer
		auto ni = clustered.neighbors.begin();
		for (auto begin = toneighs[0].begin(), end = begin; end != toneighs[0].end(); begin = end++) {
			while (*ni < begin->target) { ++ni; }
			while (end != toneighs[0].end() && begin->target == end->target) { ++end; }

			for (auto n = begin; n != end; ++n) {
				offsets.push_back(n->offset);
				coordinates.push_back(clustered.coordinates[n->offset]);
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

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, clustered.neighbors)) {
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
	searchDuplicateNodes(coordinates, [&] (esint origin, esint target) {
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
		for (size_t n = toerase.front(), last = toerase.front(), e = 0; n < clustered.noffsets.size(); ++n) {
			if (toerase.size() <= e || n != toerase[e]) {
				clustered.noffsets[last] = clustered.noffsets[n];
				clustered.coordinates[last] = clustered.coordinates[n];
				++last;
			} else {
				++e;
			}
		}
		clustered.coordinates.resize(clustered.coordinates.size() - toerase.size());
		clustered.noffsets.resize(clustered.noffsets.size() - toerase.size());
	}

	if (toinsert.size()) {
		utils::sortAndRemoveDuplicates(toinsert);
		clustered.coordinates.resize(clustered.coordinates.size() + toinsert.size());
		clustered.noffsets.resize(clustered.noffsets.size() + toinsert.size());
		for (size_t i = clustered.noffsets.size() - toinsert.size(), j = 0; j < toinsert.size(); ++i, ++j) {
			clustered.noffsets[i] = offsets[toinsert[j]];
			clustered.coordinates[i] = coordinates[toinsert[j]];
		}
	}

	searchDuplicateNodes(clustered.coordinates, [&] (esint origin, esint target) {
		clustered.nduplication.push_back({ target, origin });
	});
	std::sort(clustered.nduplication.begin(), clustered.nduplication.end());

	std::unordered_map<esint, esint> dmap;
	for (size_t n = 0; n < clustered.nduplication.size(); ++n) {
		dmap[clustered.nduplication[n].duplication] = clustered.nduplication[n].origin;
	}
}

void searchParentAndDuplicatedElements(ClusteredMesh &linked)
{
	// eoffset, etype, enodes; ....
	std::vector<std::vector<esint> > sBuffer(linked.neighbors.size()), rBuffer(linked.neighbors.size());
	if (linked.neighbors.size()) { // duplicated elements have to have all nodes held by other rank
		std::vector<int> neighbors;
		for (size_t e = 0; e < linked.eoffsets.size(); ++e) {
			neighbors.assign(linked.neighbors.begin(), linked.neighbors.end());
			for (esint n = linked.edist[e]; n < linked.edist[e + 1]; ++n) {
				size_t intersection = 0, current = 0;
				for (int r = linked.rankDistribution[linked.g2l[linked.enodes[n]]]; r < linked.rankDistribution[linked.g2l[linked.enodes[n]] + 1]; ++r) {
					while (current < neighbors.size() && neighbors[current] < linked.rankData[r]) { ++current; }
					if (current < neighbors.size() && neighbors[current] == linked.rankData[r]) {
						neighbors[intersection++] = neighbors[current++];
					}
				}
				neighbors.resize(intersection);
			}
			esint roffset = 0;
			for (auto r = neighbors.begin(); r != neighbors.end(); ++r) {
				while (linked.neighbors[roffset] < *r) { ++roffset; }
				sBuffer[roffset].push_back(linked.eoffsets[e]);
				sBuffer[roffset].push_back((esint)linked.etype[e]);
				for (esint n = linked.edist[e]; n < linked.edist[e + 1]; ++n) {
					sBuffer[roffset].push_back(linked.enodes[n]);
				}
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, linked.neighbors)) {
		eslog::internalFailure("cannot exchange duplicated elements.\n");
	}
	utils::clearVector(sBuffer);

	auto getMinimal = [&] (esint nsize, esint *nodes) {
		esint min = linked.g2l[nodes[0]];
		for (int n = 1; n < nsize; ++n) {
			esint current = linked.g2l[nodes[n]];
			if (linked.coordinates[current].x < linked.coordinates[min].x) {
				min = current;
			} else if (linked.coordinates[current].x == linked.coordinates[min].x && current < min) {
				min = current;
			}
		}
		return min;
	};

	std::vector<esint> mapDist(linked.noffsets.size() + 1);
	std::vector<esint, initless_allocator<esint> > mapData;

	for (size_t e = 0; e < linked.eoffsets.size(); ++e) {
		if (Mesh::edata[(int)linked.etype[e]].dimension == linked.dimension) {
			++mapDist[getMinimal(linked.edist[e + 1] - linked.edist[e], linked.enodes.data() + linked.edist[e])];
		} else {
			for (esint n = linked.edist[e]; n < linked.edist[e + 1]; ++n) {
				++mapDist[linked.g2l[linked.enodes[n]]];
			}
		}
	}
	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Mesh::edata[rBuffer[r][i + 1]].nodes) {
			if (Mesh::edata[rBuffer[r][i + 1]].dimension == linked.dimension) {
				++mapDist[getMinimal(Mesh::edata[rBuffer[r][i + 1]].nodes, rBuffer[r].data() + i + 2)];
			} else {
				for (int n = 0; n < Mesh::edata[rBuffer[r][i + 1]].nodes; ++n) {
					++mapDist[linked.g2l[rBuffer[r][i + 2 + n]]];
				}
			}
		}
	}
	utils::sizesToOffsets(mapDist);
	mapData.resize(mapDist.back());

	std::vector<std::pair<esint, esint> > rmap;
	std::vector<esint> _mapDist = mapDist;

	for (size_t e = 0; e < linked.eoffsets.size(); ++e) {
		if (Mesh::edata[(int)linked.etype[e]].dimension == linked.dimension) {
			mapData[_mapDist[getMinimal(linked.edist[e + 1] - linked.edist[e], linked.enodes.data() + linked.edist[e])]++] = e;
		} else {
			for (esint n = linked.edist[e]; n < linked.edist[e + 1]; ++n) {
				mapData[_mapDist[linked.g2l[linked.enodes[n]]]++] = e;
			}
		}
	}
	for (size_t r = 0; r < rBuffer.size(); ++r) {
		for (size_t i = 0; i < rBuffer[r].size(); i += 2 + Mesh::edata[rBuffer[r][i + 1]].nodes) {
			if (Mesh::edata[rBuffer[r][i + 1]].dimension == linked.dimension) {
				mapData[_mapDist[getMinimal(Mesh::edata[rBuffer[r][i + 1]].nodes, rBuffer[r].data() + i + 2)]++] = linked.eoffsets.size() + rmap.size();
				rmap.push_back({ (esint)r, (esint)i });
			} else {
				for (int n = 0; n < Mesh::edata[rBuffer[r][i + 1]].nodes; ++n) {
					mapData[_mapDist[linked.g2l[rBuffer[r][i + 2 + n]]]++] = linked.eoffsets.size() + rmap.size();
				}
				rmap.push_back({ (esint)r, (esint)i });
			}
		}
	}
	utils::clearVector(_mapDist);

	struct __parent__ {
		esint boundary;
		esint parent;
	};

	std::vector<__parent__> parents;
	std::vector<esint, initless_allocator<esint> > shared, tmpShared;
	auto merge = [&] (esint nsize, esint *nodes) {
		esint min = linked.g2l[nodes[0]];
		shared.assign(mapData.begin() + mapDist[min], mapData.begin() + mapDist[min + 1]);
		for (esint n = 1; n < nsize; ++n) {
			esint current = linked.g2l[nodes[n]];
			tmpShared.swap(shared);
			shared.resize(tmpShared.size() + mapDist[current + 1] - mapDist[current]);
			std::merge(mapData.begin() + mapDist[current], mapData.begin() + mapDist[current + 1], tmpShared.begin(), tmpShared.end(), shared.begin());
			if (linked.coordinates[current].x < linked.coordinates[min].x) {
				min = current;
			} else if (linked.coordinates[current].x == linked.coordinates[min].x && current < min) {
				min = current;
			}
		}
		return min;
	};

	std::unordered_set<esint> duplicated;
	std::vector<esint> _checkBuffer;
	_checkBuffer.reserve(2 * 20);
	auto areSame = [&] (esint size1, esint *nodes1, esint size2, esint *nodes2) {
		if (size1 == size2) {
			_checkBuffer.clear();
			for (esint n = 0; n < size1; ++n) { _checkBuffer.push_back(linked.g2l[nodes1[n]]); }
			for (esint n = 0; n < size2; ++n) { _checkBuffer.push_back(linked.g2l[nodes2[n]]); }
			std::sort(_checkBuffer.begin(), _checkBuffer.begin() + size1);
			std::sort(_checkBuffer.begin() + size1, _checkBuffer.end());
			return memcmp(_checkBuffer.data(), _checkBuffer.data() + size1, sizeof(esint) * size1) == 0;
		}
		return false;
	};
	for (size_t e1 = 0; e1 < linked.eoffsets.size(); ++e1) {
		if (Mesh::edata[(int)linked.etype[e1]].dimension == linked.dimension) {
			esint min = merge(linked.edist[e1 + 1] - linked.edist[e1], linked.enodes.data() + linked.edist[e1]);
			if (duplicated.count(e1) == 0) {
				for (esint n = mapDist[min]; n < mapDist[min + 1]; ++n) {
					size_t e2 = mapData[n];
					if (e1 < e2) { // if e2 is lower, the pair was already processed
						if (e2 < linked.eoffsets.size() && Mesh::element(linked.etype[e2]).dimension == linked.dimension) { // local
							if (areSame(linked.edist[e1 + 1] - linked.edist[e1], linked.enodes.data() + linked.edist[e1], linked.edist[e2 + 1] - linked.edist[e2], linked.enodes.data() + linked.edist[e2])) {
								linked.eduplication.push_back({ linked.eoffsets[e1], linked.eoffsets[e2] });
								duplicated.insert(e2);
							}
						} else { // from neighbors
							const auto &r = rmap[e2 - linked.eoffsets.size()];
							if (Mesh::element(rBuffer[r.first][r.second + 1]).dimension == linked.dimension) {
								if (areSame(linked.edist[e1 + 1] - linked.edist[e1], linked.enodes.data() + linked.edist[e1], Mesh::element(rBuffer[r.first][r.second + 1]).nodes, rBuffer[r.first].data() + r.second + 2)) {
									linked.eduplication.push_back({ linked.eoffsets[e1], rBuffer[r.first][r.second] });
									duplicated.insert(e2);
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
					if (e2 < linked.eoffsets.size()) {
						if (Mesh::element(linked.etype[e2]).nodes == end - begin) {
							parents.push_back({ linked.eoffsets[e2], linked.eoffsets[e1] });
						}
					} else {
						auto &roffset = rmap[e2 - linked.eoffsets.size()];
						if (Mesh::element(rBuffer[roffset.first][roffset.second + 1]).nodes == end - begin) {
							parents.push_back({ rBuffer[roffset.first][roffset.second], linked.eoffsets[e1] });
						}
					}
				}
			}
		}
	}
	std::sort(linked.eduplication.begin(), linked.eduplication.end());
	// process parents
}

}
}
