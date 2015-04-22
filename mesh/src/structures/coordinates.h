#ifndef COORDINATES_H_
#define COORDINATES_H_

#include "../elements/1D/point.h"
#include "../loader.h"

#include <vector>

class Coordinates
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Coordinates &c);

	Coordinates(): _points(0), _offset(0) { };
	Coordinates(const char *fileName);
	Coordinates(size_t size, idx_t offset): _points(size), _offset(offset) {};

	void resize(size_t size)
	{
		_points.resize(size);
	}

	/** @brief Index correction. C/C++ indexes start at 0 while first Point is indexed by 1. */
	void setOffset(idx_t offset)
	{
		_offset = offset;
	}

	idx_t getOffset() const
	{
		return _offset;
	}

	size_t size() const
	{
		return _points.size();
	}

	const Point& operator[](idx_t index) const
	{
		return _points[index - _offset];
	}

	Point& operator[](idx_t index)
	{
		return _points[index - _offset];
	}

private:
	std::vector<Point> _points;

	/** @brief Correction between C/C++ and Point indexing. */
	idx_t _offset;
};



#endif /* COORDINATES_H_ */
