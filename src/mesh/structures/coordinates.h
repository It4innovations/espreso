#ifndef COORDINATES_H_
#define COORDINATES_H_

#include "esbasis.h"
#include "../settings/setting.h"

#include <algorithm>
#include <vector>
#include <map>
#include <fstream>

namespace espreso
{

class Coordinates
{
	friend class Mesh;
	friend std::ostream& operator<<(std::ostream& os, const Coordinates &c);

public:
	Coordinates(): _points(0), _clusterIndex(1) { };

	void add(const Point &point, eslocal clusterIndex, esglobal globalIndex)
	{
		_points.push_back(point);
		_clusterIndex[0].push_back(clusterIndex);
		_globalIndex.push_back(globalIndex);
		_globalMap[globalIndex] = clusterIndex;
	}

	void reserve(size_t size)
	{
		_points.reserve(size);
		_globalIndex.reserve(size);
		_clusterIndex[0].reserve(size);
	}

	void clear()
	{
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
		return _globalMap.find(index)->second;
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
		return _points.size();
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

private:
	std::vector<Point> _points;

	/** @brief Local point to cluster index. */
	std::vector<std::vector<eslocal> > _clusterIndex;

	/** @brief Point to global index */
	std::vector<esglobal> _globalIndex;
	std::map<esglobal, eslocal> _globalMap;
};

}

#endif /* COORDINATES_H_ */
