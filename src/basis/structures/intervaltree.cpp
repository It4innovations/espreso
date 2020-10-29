
#include "intervaltree.h"
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdio>

using namespace espreso;

IntervalTree::IntervalTree(const std::vector<Point> &start, const std::vector<Point> &end)
: istart(start), iend(end)
{
	if (end.size()) {
		min = start.front();
		max = end.front();
		for (size_t i = 1; i < end.size(); ++i) {
			min.x = std::min(start[i].x, min.x);
			min.y = std::min(start[i].y, min.y);
			min.z = std::min(start[i].z, min.z);
			max.x = std::max(end[i].x, max.x);
			max.y = std::max(end[i].y, max.y);
			max.z = std::max(end[i].z, max.z);
		}
	}
	build();
}

void IntervalTree::build()
{
	struct comp {
		comp(const std::vector<Point> &vals, int d): vals(vals), d(d) {}
		bool operator()(esint i, esint j) { return vals[i][d] < vals[j][d]; }

		const std::vector<Point> &vals;
		int d;
	};

	size_t bucketsize = 8;
	size = max - min;
	permutation.resize(iend.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	levels = iend.size() < bucketsize ? 0 : std::floor(std::log2(iend.size() / bucketsize));
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
				return iend[i][splitters[index].d] < iend[j][splitters[index].d];
			});
			if (end - begin) {
				splitters[index].end = iend[permutation[splitters[index].index]][splitters[index].d];
			}
			while ( // move to the last coordinate with the same value as mid
					splitters[index].index + 1 < end &&
					splitters[index].end == iend[permutation[splitters[index].index + 1]][splitters[index].d]) {
				++splitters[index].index;
			}
			for (esint c = splitters[index].index + 2; c < end; ++c) { // there can be another in the rest of array
				if (splitters[index].end == iend[permutation[c]][splitters[index].d]) {
					std::swap(permutation[++splitters[index].index], permutation[c--]);
				}
			}
			splitters[index].index = std::min(splitters[index].index + 1, end);

			splitters[index].start = splitters[index].end;
			for (esint c = splitters[index].index; c < end; ++c) {
				splitters[index].start = std::min(splitters[index].start, istart[permutation[c]][splitters[index].d]);
			}
		}
	}
}
