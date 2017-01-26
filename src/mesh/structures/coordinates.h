#ifndef COORDINATES_H_
#define COORDINATES_H_

#include <vector>
#include <algorithm>

#include "../../basis/point/point.h"

namespace espreso
{

struct G2L {
	esglobal global;
	esglobal local;

	G2L(esglobal global, esglobal local): global(global), local(local) {};

	bool operator<(const G2L &other) const
	{
		return global < other.global;
	}
};

struct Coordinates
{
	friend class Mesh;
	friend std::ostream& operator<<(std::ostream& os, const Coordinates &c);

public:
	Coordinates(): _clusterSize(0), _clusterIndex(1) { };

	void resize(size_t size)
	{
		_clusterSize = size;
		_points.reserve(size);
		_clusterIndex.resize(1);
		_clusterIndex[0].resize(size);
	}

	void add(const Point &point, eslocal clusterIndex, esglobal globalIndex)
	{
		_points.push_back(point);
		_clusterIndex[0].push_back(clusterIndex);
		_globalIndex.push_back(globalIndex);
		_globalMapping.push_back(G2L(globalIndex, clusterIndex));
	}

	void reserve(size_t size)
	{
		_points.reserve(size);
		_globalIndex.reserve(size);
		_clusterIndex[0].reserve(size);
	}

	void clear()
	{
		_clusterSize = 0;
		_points.clear();
		_globalIndex.clear();
		_clusterIndex.resize(1);
		_clusterIndex[0].clear();
	}

	const Point& get(eslocal index, eslocal part) const
	{
		return _points[_clusterIndex[part][index]];
	}

	eslocal clusterIndex(eslocal index, eslocal part) const
	{
		return _clusterIndex[part][index];
	}

	esglobal globalIndex(eslocal index, eslocal part) const
	{
		return _globalIndex[_clusterIndex[part][index]];
	}

	esglobal globalIndex(eslocal index) const
	{
		return _globalIndex[index];
	}

	eslocal clusterIndex(esglobal index) const
	{
		return std::lower_bound(_globalMapping.begin(), _globalMapping.end(), index, [] (const G2L &mapping, esglobal index) {
			return mapping.global < index;
		})->local;
	}

	eslocal localIndex(eslocal index, eslocal part) const
	{
		std::vector< eslocal >::const_iterator it;
		it = lower_bound(_clusterIndex[part].begin(), _clusterIndex[part].end(), index) ;

		if (it == _clusterIndex[part].end() || *it != index) {
			return -1;
		} else {
			return it - _clusterIndex[part].begin();
		}
	}


	size_t clusterSize() const
	{
		return _clusterSize < _points.size() ? _points.size() : _clusterSize;
	}

	size_t localSize(eslocal part) const
	{
		return _clusterIndex[part].size();
	}

	size_t parts() const
	{
		return _clusterIndex.size();
	}

	const std::vector<esglobal>& clusterToGlobal() const
	{
		return _globalIndex;
	}

	const std::vector<eslocal>& localToCluster(eslocal part) const
	{
		return _clusterIndex[part];
	}

	void localClear()
	{
		_clusterIndex.clear();
	}

	void localResize(size_t size)
	{
		_clusterIndex.resize(size);
	}

	const Point& operator[](eslocal index) const
	{
		return _points[index];
	}

	Point& operator[](eslocal index)
	{
		return _points[index];
	}

// private:
	size_t _clusterSize;
	std::vector<Point> _points;

	/** @brief Local point to cluster index. */
	std::vector<std::vector<eslocal> > _clusterIndex;

	/** @brief Point to global index */
	std::vector<esglobal> _globalIndex;
	std::vector<G2L> _globalMapping;
};

}

#endif /* COORDINATES_H_ */
