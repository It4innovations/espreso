#ifndef COORDINATES_H_
#define COORDINATES_H_

#include "../elements/1D/point.h"
#include "../loader.h"

#include <vector>

namespace mesh {

class Coordinates
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Coordinates &c);

	Coordinates(): _points(0), _clusterIndex(1) { };
	Coordinates(const char *fileName);

	void add(const Point &point, esint clusterIndex, eslong globalIndex)
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

	const Point& get(esint index, esint part) const
	{
		return _points[_clusterIndex[part][index]];
	}

	esint clusterIndex(esint index, esint part) const
	{
		return _clusterIndex[part][index];
	}

	eslong globalIndex(esint index, esint part) const
	{
		return _globalIndex[_clusterIndex[part][index]];
	}

	esint clusterSize() const
	{
		return _points.size();
	}

	esint localSize(esint part) const
	{
		return _clusterIndex[part].size();
	}

	size_t size() const
	{
		return _points.size();
	}

	const std::vector<esint>& localToCluster(esint part) const
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

	void computeLocal(esint part, std::vector<esint> &nodeMap, size_t size);

	double* data()
	{
		void *tmp = &_points[0];
		return static_cast<double*>(tmp);
	}

	const Point& operator[](esint index) const
	{
		return _points[index];
	}

	Point& operator[](esint index)
	{
		return _points[index];
	}


private:
	std::vector<Point> _points;

	/** @brief Local point to cluster index. */
	std::vector<std::vector<esint> > _clusterIndex;

	/** @brief Point to global index */
	std::vector<eslong> _globalIndex;
};

}

#endif /* COORDINATES_H_ */
