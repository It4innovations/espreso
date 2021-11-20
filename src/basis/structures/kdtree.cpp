
#include "kdtree.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdio>
#include <functional>
#include <utility>

namespace espreso {

template <typename T>
static void _build(KDTree<T> *tree, size_t bucketsize)
{
	auto compxyz = [&] (esint i, esint j) {
		if (tree->coordinates[i].x == tree->coordinates[j].x) {
			if (tree->coordinates[i].y == tree->coordinates[j].y) {
				return tree->coordinates[i].z < tree->coordinates[j].z;
			}
			return tree->coordinates[i].y < tree->coordinates[j].y;
		}
		return tree->coordinates[i].x < tree->coordinates[j].x;
	};

	tree->size = tree->max - tree->min;
	tree->permutation.resize(tree->coordinates.size());
	std::iota(tree->permutation.begin(), tree->permutation.end(), 0);
	tree->levels = tree->coordinates.size() < bucketsize ? 0 : std::floor(std::log2(tree->coordinates.size() / bucketsize));
	tree->splitters.resize(std::exp2(tree->levels));

	_Point<T> box = tree->size; // uniform division (denoted by the level)
	for (esint ll = 0, intervals = 1; ll < tree->levels; ++ll, intervals *= 2) {
		int splitter = box.x < box.y ? box.y < box.z ? 2 : 1 : box.x < box.z ? 2 : 0;
		box[splitter] /= 2;

		if(ll + 1 == tree->levels){
			tree->leaf_intervals.resize(intervals);
		}

		#pragma omp parallel for
		for (esint i = 0; i < intervals; ++i) {
			esint index = std::exp2(ll) + i;
			esint begin = tree->begin(index);
			esint end = tree->end(index);

			if(ll + 1 == tree->levels){
				tree->leaf_intervals[i].begin = begin;
				tree->leaf_intervals[i].end = end;
			}

			tree->splitters[index].d = splitter;
			tree->splitters[index].index = begin + (end - begin) / 2;
			

			std::nth_element(tree->permutation.begin() + begin, tree->permutation.begin() + tree->splitters[index].index, tree->permutation.begin() + end, [&] (esint i, esint j) {
				return tree->coordinates[i][tree->splitters[index].d] < tree->coordinates[j][tree->splitters[index].d];
			});
			if (end - begin) {
				tree->splitters[index].value = tree->coordinates[tree->permutation[tree->splitters[index].index]][tree->splitters[index].d];
			}
			while ( // move to the last coordinate with the same value as mid
					tree->splitters[index].index + 1 < end &&
					tree->coordinates[tree->permutation[tree->splitters[index].index]][tree->splitters[index].d] == tree->coordinates[tree->permutation[tree->splitters[index].index + 1]][tree->splitters[index].d]) {
				++tree->splitters[index].index;
			}
			for (esint c = tree->splitters[index].index + 2; c < end; ++c) { // there can be another in the rest of array
				if (tree->coordinates[tree->permutation[tree->splitters[index].index]][tree->splitters[index].d] == tree->coordinates[tree->permutation[c]][tree->splitters[index].d]) {
					std::swap(tree->permutation[++tree->splitters[index].index], tree->permutation[c--]);
				}
			}
			tree->splitters[index].index = std::min(tree->splitters[index].index + 1, end);
			if (ll + 1 == tree->levels) {
				std::sort(tree->permutation.begin() + begin, tree->permutation.begin() + tree->splitters[index].index, compxyz);
				std::sort(tree->permutation.begin() + tree->splitters[index].index, tree->permutation.begin() + end, compxyz);
			}
		}
	}
	if (tree->levels == 0) {
		std::sort(tree->permutation.begin(), tree->permutation.end(), compxyz);
	}
}

template <> void KDTree<float>::build(size_t bucketsize) { _build(this, bucketsize); }
template <> void KDTree<double>::build(size_t bucketsize) { _build(this, bucketsize); }

}
