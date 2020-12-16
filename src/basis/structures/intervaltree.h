
#ifndef SRC_BASIS_STRUCTURES_INTERVALTREE_H_
#define SRC_BASIS_STRUCTURES_INTERVALTREE_H_

#include "basis/containers/point.h"

#include <vector>

namespace espreso {

class IntervalTree {
	struct splitter {
		int d;
		esint index;
		double start, end;
	};

public:
	IntervalTree(const std::vector<Point> &start, const std::vector<Point> &end);

	esint begin(esint interval)
	{
		esint left = interval;
		while (left % 2 == 0) { left /= 2; }
		return left == 1 ? 0 : splitters[left / 2].index;
	}

	esint end(esint interval)
	{
		esint right = interval;
		while (right % 2 == 1) { right /= 2; }
		return right == 0 ? permutation.size() : splitters[right / 2].index;
	}

	void boxMin(esint interval, Point &min)
	{
		min = this->min;
		esint parent = interval / 2;
		while (parent) {
			if (interval % 2 == 1) {
				if (min[splitters[parent].d] < splitters[parent].end) {
					min[splitters[parent].d] = splitters[parent].end;
				}
			}
			interval = parent;
			parent /= 2;
		}
	}

	void traverse(esint node, const Point &start, const Point &end, std::vector<esint> &intervals)
	{
		if (node < (esint)splitters.size()) {
			if (start[splitters[node].d] < splitters[node].end) {
				traverse(2 * node, start, end, intervals);
			}
			if (splitters[node].start < end[splitters[node].d]) {
				traverse(2 * node + 1, start, end, intervals);
			}
		} else {
			intervals.push_back(node);
		}
	}

	const std::vector<Point> &istart, &iend;

	esint levels;
	Point min, max, size;
	std::vector<esint> permutation;
	std::vector<splitter> splitters;

protected:
	void build();
};

}

#endif /* SRC_BASIS_STRUCTURES_INTERVALTREE_H_ */
