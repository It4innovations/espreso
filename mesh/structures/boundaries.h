#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include <sstream>

#include "../elements/elements.h"
#include "cilk/cilk.h"

namespace mesh {

class Boundaries
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Boundaries &f);

	void clear()
	{
		_boundaries.clear();
		_corners.clear();
	}

	void resize(size_t size)
	{
		_boundaries.resize(size);
		_corners.resize(size, false);
	}

	size_t size() const
	{
		return _boundaries.size();
	}

	const std::set<eslocal>& operator[](size_t position) const
	{
		return _boundaries[position];
	}

	std::set<eslocal>& operator[](size_t position)
	{
		return _boundaries[position];
	}

	void setCorner(size_t index)
	{
		_corners[index] = true;
	}

	bool isCorner(size_t index) const
	{
		return _corners[index];
	}

	// prepare for future improvements
	eslocal index(size_t position) const {
		return position;
	}

private:
	/** @brief Keeps mapping of nodes to mesh parts. */
	std::vector<std::set<eslocal> > _boundaries;

	/** @brief Keeps information whether a point is the corner. */
	std::vector<bool> _corners;
};

}

#endif /* BOUNDARIES_H_ */
