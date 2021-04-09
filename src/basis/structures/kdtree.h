
#ifndef SRC_BASIS_STRUCTURES_KDTREE_H_
#define SRC_BASIS_STRUCTURES_KDTREE_H_

#include "basis/containers/point.h"

#include <vector>
#include <cstdio>

namespace espreso {

class KDTree {
	enum: size_t { defaultBucketSize = 8 };

	struct splitter {
		int d;
		esint index;
		double value;
	};

	struct leaf_intervals {
		esint begin;
		esint end;
	};

public:
	KDTree(const std::vector<Point> &coordinates, size_t bucketsize = defaultBucketSize);
	KDTree(const std::vector<Point> &coordinates, Point &min, Point &max);

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

	esint getMidPoint(esint interval)
	{
		esint i;
		double dx, dy, dz, d, dmin = 0.0f;

		esint b = this->leaf_intervals[interval].begin;
		esint e = this->leaf_intervals[interval].end;

		esint out = -1;

		Point centroid;

		for(i = b; i < e; ++i){
			centroid += this->coordinates[this->permutation[i]];
		}
		centroid /= (e-b);

		for(i = b; i < e; ++i){
			dx = this->coordinates[this->permutation[i]][0] - centroid[0];
			dy = this->coordinates[this->permutation[i]][1] - centroid[1];
			dz = this->coordinates[this->permutation[i]][2] - centroid[2];
			d = dx*dx + dy*dy + dz*dz;

			if(d < dmin || out == -1){
				dmin = d;
				out = this->permutation[i];
			}
		}

		// printf("Bucket #%3d of size %d, Centroid: (%f, %f, %f)\n", interval, e - b, centroid.x, centroid.y, centroid.z);
		// printf("The most central bucket member: #%d (%f, %f, %f)\n", out, this->coordinates[out].x, this->coordinates[out].y, this->coordinates[out].z);

		return out;
	};

	esint getNLeaves(){
		return this->leaf_intervals.size();
	};

	void boxMin(esint interval, Point &min) {
		min = this->min;
		esint parent = interval / 2;
		while (parent) {
			if (interval % 2 == 1) {
				if (min[splitters[parent].d] < splitters[parent].value) {
					min[splitters[parent].d] = splitters[parent].value;
				}
			}
			interval = parent;
			parent /= 2;
		}
	}

	const std::vector<Point> &coordinates;

	esint levels;
	Point min, max, size;
	std::vector<esint> permutation;
	std::vector<splitter> splitters;
	std::vector<leaf_intervals> leaf_intervals;

protected:
	void build(size_t bucketsize = defaultBucketSize);
};

}

#endif /* SRC_BASIS_STRUCTURES_KDTREE_H_ */
