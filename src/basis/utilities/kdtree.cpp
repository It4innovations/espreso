
#include "kdtree.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdio>
#include <functional>
#include <utility>

using namespace espreso;

KDTree::KDTree(const std::vector<Point> &coordinates)
: coordinates(coordinates)
{
	if (coordinates.size()) {
		min = coordinates.front();
		max = coordinates.front();
		for (size_t i = 1; i < coordinates.size(); ++i) {
			min.x = std::min(coordinates[i].x, min.x);
			min.y = std::min(coordinates[i].y, min.y);
			min.z = std::min(coordinates[i].z, min.z);
			max.x = std::max(coordinates[i].x, max.x);
			max.y = std::max(coordinates[i].y, max.y);
			max.z = std::max(coordinates[i].z, max.z);
		}
	}
	build();
}

KDTree::KDTree(const std::vector<Point> &coordinates, Point &min, Point &max)
: coordinates(coordinates), min(min), max(max)
{
	build();
}

void KDTree::build()
{
	size_t bucketsize = 8;
	auto compxyz = [&] (esint i, esint j) {
		if (coordinates[i].x == coordinates[j].x) {
			if (coordinates[i].y == coordinates[j].y) {
				return coordinates[i].z < coordinates[j].z;
			}
			return coordinates[i].y < coordinates[j].y;
		}
		return coordinates[i].x < coordinates[j].x;
	};

	size = max - min;
	permutation.resize(coordinates.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	levels = coordinates.size() < bucketsize ? 0 : std::floor(std::log2(coordinates.size() / bucketsize));
	splitters.resize(std::exp2(levels));

	Point box = size; // uniform division (denoted by the level)
	for (esint ll = 0, intervals = 1; ll < levels; ++ll, intervals *= 2) {
		int splitter = box.x < box.y ? box.y < box.z ? 2 : 1 : box.x < box.z ? 2 : 0;
		box[splitter] /= 2;

		#pragma omp parallel for
		for (esint i = 0; i < intervals; ++i) {
			esint index = std::exp2(ll) + i;
			esint begin = this->begin(index);
			esint end = this->end(index);

			splitters[index].d = splitter;
			splitters[index].index = begin + (end - begin) / 2;
			std::nth_element(permutation.begin() + begin, permutation.begin() + splitters[index].index, permutation.begin() + end, [&] (esint i, esint j) {
				return coordinates[i][splitters[index].d] < coordinates[j][splitters[index].d];
			});
			if (end - begin) {
				splitters[index].value = coordinates[permutation[splitters[index].index]][splitters[index].d];
			}
			while ( // move to the last coordinate with the same value as mid
					splitters[index].index + 1 < end &&
					coordinates[permutation[splitters[index].index]][splitters[index].d] == coordinates[permutation[splitters[index].index + 1]][splitters[index].d]) {
				++splitters[index].index;
			}
			for (esint c = splitters[index].index + 2; c < end; ++c) { // there can be another in the rest of array
				if (coordinates[permutation[splitters[index].index]][splitters[index].d] == coordinates[permutation[c]][splitters[index].d]) {
					std::swap(permutation[++splitters[index].index], permutation[c--]);
				}
			}
			splitters[index].index = std::min(splitters[index].index + 1, end);
			if (ll + 1 == levels) {
				std::sort(permutation.begin() + begin, permutation.begin() + splitters[index].index, compxyz);
				std::sort(permutation.begin() + splitters[index].index, permutation.begin() + end, compxyz);
			}
		}
	}
	if (levels == 0) {
		std::sort(permutation.begin(), permutation.end(), compxyz);
	}
}
