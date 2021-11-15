
#ifndef SRC_BASIS_SFC_SPACEFILLINGCURVE_H_
#define SRC_BASIS_SFC_SPACEFILLINGCURVE_H_

#include "basis/containers/point.h"

#include <cstddef>
#include <cmath>
#include <vector>
#include <functional>

namespace espreso {

// if SFCDEPTH > 10 => set buckets to size_t
#define SFCDEPTH 10

struct SpaceFillingCurve {

	SpaceFillingCurve(size_t dimension, size_t depth, size_t npoints, Point* coordinates);
	virtual ~SpaceFillingCurve() {}

	size_t depth() const { return _depth; }
	size_t dimension() const { return _dimension; }
	size_t buckets(size_t depth) const { return std::pow((size_t)1 << depth, _dimension); }
	size_t bucketSize() const { return _dimension == 2 ? 4 : 8; }

	// computed from splitters
	void addSFCNeighbors(size_t depth, size_t index, std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors);
	void addXYNeighbors(size_t depth, size_t x, size_t y, std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors);
	void addXYZNeighbors(size_t depth, size_t x, size_t y, size_t z, std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors);

	size_t getBucket(const Point &p) const { return _dimension == 2 ? D2toD1(p) : D3toD1(p); }
	size_t getBucket(const double &value, int d) { return std::floor(_n * (value - _origin[d]) / _size[d]);}

protected:
	size_t D2toD1(const Point &p) const
	{
		size_t x = std::floor(_n * (p.x - _origin.x) / _size.x);
		size_t y = std::floor(_n * (p.y - _origin.y) / _size.y);

		return D2toD1(_n, x, y);
	}
	size_t D3toD1(const Point &p) const
	{
		size_t x = std::floor(_n * (p.x - _origin.x) / _size.x);
		size_t y = std::floor(_n * (p.y - _origin.y) / _size.y);
		size_t z = std::floor(_n * (p.z - _origin.z) / _size.z);

		return D3toD1(_n, x, y, z);
	}

	virtual void D1toD2(size_t n, size_t d, size_t &x, size_t &y) const = 0;
	virtual void D1toD3(size_t n, size_t d, size_t &x, size_t &y, size_t &z) const = 0;
	virtual size_t D2toD1(size_t n, size_t x, size_t y) const = 0;
	virtual size_t D3toD1(size_t n, size_t x, size_t y, size_t z) const = 0;

	size_t _dimension;
	size_t _depth, _n;
	Point _origin, _size;
private:
	std::pair<size_t, size_t> getXYZBucket(size_t depth, size_t x, size_t y, size_t z);
};

}


#endif /* SRC_BASIS_SFC_SPACEFILLINGCURVE_H_ */
