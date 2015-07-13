#ifndef COORDINATES_H_
#define COORDINATES_H_

#include "../elements/1D/point.h"
#include "../loader.h"

#include <algorithm>
#include <vector>
#include <map>

namespace mesh
{

class CoordinatesProperty
{
public:
	friend std::ostream& operator<<(std::ostream& os, const CoordinatesProperty &cp);

	CoordinatesProperty() { };
	void load(const char *fileName);

	double& operator[](eslocal index)
	{
		return _mapping[index];
	}

	const std::map<eslocal, double>& values() const
	{
		return _mapping;
	}

private:
	std::map<eslocal, double> _mapping;
};

class Coordinates
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Coordinates &c);

	Coordinates(): _points(0), _clusterIndex(1), _property(CP::SIZE) { };
	Coordinates(const char *fileName);
	Coordinates(const Ansys &a);

	void add(const Point &point, eslocal clusterIndex, esglobal globalIndex)
	{
		_points.push_back(point);
		_clusterIndex[0].push_back(clusterIndex);
		_globalIndex.push_back(globalIndex);
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
		return lower_bound(_globalIndex.begin(), _globalIndex.end(), index) - _globalIndex.begin();
	}

	eslocal localIndex(eslocal index, eslocal part) const
	{

		std::vector< eslocal >::const_iterator it;
		it = lower_bound(_clusterIndex[part].begin(), _clusterIndex[part].end(), index) ;

		if (it == _clusterIndex[part].end() || *it != index)
			return -1;
		else
			return it - _clusterIndex[part].begin();
	}


	eslocal clusterSize() const
	{
		return _points.size();
	}

	eslocal localSize(eslocal part) const
	{
		return _clusterIndex[part].size();
	}

	size_t size() const
	{
		return _points.size();
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

	void computeLocal(eslocal part, std::vector<eslocal> &nodeMap, size_t size);

	double* data()
	{
		void *tmp = &_points[0];
		return static_cast<double*>(tmp);
	}

	const Point& operator[](eslocal index) const
	{
		return _points[index];
	}

	Point& operator[](eslocal index)
	{
		return _points[index];
	}

	CoordinatesProperty& property(CP::Property property)
	{
		return _property[property];
	}

	const CoordinatesProperty& property(CP::Property property) const
	{
		return _property[property];
	}

private:
	void readFromFile(const char *fileName);

	std::vector<Point> _points;

	/** @brief Local point to cluster index. */
	std::vector<std::vector<eslocal> > _clusterIndex;

	/** @brief Point to global index */
	std::vector<esglobal> _globalIndex;

	/** @brief Coordinates properties */
	std::vector<CoordinatesProperty> _property;
};

}

#endif /* COORDINATES_H_ */
