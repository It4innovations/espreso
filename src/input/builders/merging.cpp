
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

namespace espreso {
namespace builder {

static void searchDuplicateNodes(std::vector<_Point<esfloat> > &coordinates, std::function<void(esint origin, esint target)> merge)
{
	KDTree<esfloat> tree(coordinates);

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

void mergeDuplicatedNodes(const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered)
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
	std::vector<_Point<esfloat> > coordinates;
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

//	std::vector<std::pair<esint, esint> > buckets;
	if (toinsert.size()) {
		utils::sortAndRemoveDuplicates(toinsert);
		clustered.coordinates.resize(clustered.coordinates.size() + toinsert.size());
		clustered.noffsets.resize(clustered.noffsets.size() + toinsert.size());
		for (size_t i = clustered.noffsets.size() - toinsert.size(), j = 0; j < toinsert.size(); ++i, ++j) {
			clustered.noffsets[i] = offsets[toinsert[j]];
			clustered.coordinates[i] = coordinates[toinsert[j]];
//			buckets.push_back(std::make_pair(clustered.noffsets[i], sfc.getBucket(clustered.coordinates[i])));
//			for (size_t r = 0; r < _nregsize; ++r) {
//				_nregions[i * _nregsize + r] = regions[toinsert[j] * _nregsize + r];
//			}
		}
	}
//	{
//		// update original buckets
//		Communication::allGatherUnknownSize(buckets);
//		for (size_t i = 0; i < buckets.size(); ++i) {
//			if (_nDistribution[info::mpi::rank] <= buckets[i].first && buckets[i].first < _nDistribution[info::mpi::rank + 1]) {
//				auto it = std::lower_bound(_nIDs.begin(), _nIDs.end(), buckets[i].first);
//				if (it != _nIDs.end() && *it == buckets[i].first) {
//					_nBuckets[it - _nIDs.begin()] = buckets[i].second;
//				} else {
//					eslog::internalFailure("cannot update duplicated node bucket.\n");
//				}
//			}
//		}
//	}
//
	// we need to sum regions for all duplicate nodes
	// there can be more optimal solution, but this keep the rest of the code unchanged
//	searchDuplicateNodes(clustered.coordinates, [&] (esint origin, esint target) {
////		clustered._duplicateNodes.push_back(MeshBuilder::Duplicate{ clustered.noffsets[origin], clustered.noffsets[target], origin, target });
//	});
//	std::sort(clustered._duplicateNodes.begin(), clustered._duplicateNodes.end(), MeshBuilder::Duplicate());
//	for (auto dup = _meshData._duplicateNodes.begin(); dup != _meshData._duplicateNodes.end(); ++dup) {
//		for (size_t i = 0; i < _nregsize; i++) {
//			_nregions[_nregsize * dup->toffset + i] |= _nregions[_nregsize * dup->idoffset + i];
//		}
//	}
//	for (auto dup = _meshData._duplicateNodes.begin(); dup != _meshData._duplicateNodes.end(); ++dup) {
//		for (size_t i = 0; i < _nregsize; i++) {
//			_nregions[_nregsize * dup->idoffset + i] |= _nregions[_nregsize * dup->toffset + i];
//		}
//	}
}

}
}
