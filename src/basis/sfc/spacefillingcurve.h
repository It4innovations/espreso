
#ifndef SRC_BASIS_SFC_SPACEFILLINGCURVE_H_
#define SRC_BASIS_SFC_SPACEFILLINGCURVE_H_

#include "basis/containers/allocators.h"
#include "basis/containers/point.h"

#include <cstddef>
#include <cmath>
#include <vector>
#include <functional>

namespace espreso {

// if SFCDEPTH > 10 => set buckets to size_t
#define SFCDEPTH 10

template <typename T>
struct SpaceFillingCurve {

	SpaceFillingCurve(size_t dimension, size_t depth, size_t npoints, _Point<T>* coordinates);
	virtual ~SpaceFillingCurve() {}

	size_t buckets(size_t depth) const { return std::pow((size_t)1 << depth, dimension); }
	size_t bucketSize() const { return dimension == 2 ? 4 : 8; }

	// computed from splitters
	void addSFCNeighbors(size_t depth, size_t index, const ivector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const;
	void addXYNeighbors(size_t depth, size_t x, size_t y, const ivector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const;
	void addXYZNeighbors(size_t depth, size_t x, size_t y, size_t z, const ivector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const;

	size_t getBucket(const _Point<T> &p) const { return dimension == 2 ? D2toD1(p) : D3toD1(p); }
	size_t getBucket(const double &value, int d) const { return std::floor(n * (value - origin[d]) / size[d]);}

	size_t D2toD1(const _Point<T> &p) const
	{
		size_t x = std::floor(n * (p.x - origin.x) / size.x);
		size_t y = std::floor(n * (p.y - origin.y) / size.y);

		return D2toD1(n, x, y);
	}
	size_t D3toD1(const _Point<T> &p) const
	{
		size_t x = std::floor(n * (p.x - origin.x) / size.x);
		size_t y = std::floor(n * (p.y - origin.y) / size.y);
		size_t z = std::floor(n * (p.z - origin.z) / size.z);

		return D3toD1(n, x, y, z);
	}

	virtual void D1toD2(size_t n, size_t d, size_t &x, size_t &y) const = 0;
	virtual void D1toD3(size_t n, size_t d, size_t &x, size_t &y, size_t &z) const = 0;
	virtual size_t D2toD1(size_t n, size_t x, size_t y) const = 0;
	virtual size_t D3toD1(size_t n, size_t x, size_t y, size_t z) const = 0;

	size_t dimension;
	size_t depth, n;
	_Point<T> origin, size;
};

}


#endif /* SRC_BASIS_SFC_SPACEFILLINGCURVE_H_ */
