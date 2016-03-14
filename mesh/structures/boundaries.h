#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include <sstream>

#include "../elements/elements.h"

namespace espreso {

class Boundaries
{

public:

	friend class Mesh;
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

	const std::vector<eslocal>& operator[](size_t position) const
	{
		return _boundaries[position];
	}

	std::vector<eslocal>& operator[](size_t position)
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

	bool isAveraging(eslocal corner) const
	{
		return _averaging.find(corner) != _averaging.end();
	}

	const std::vector<eslocal>& averaging(eslocal corner) const
	{
		return _averaging.find(corner)->second;
	}

	std::vector<eslocal>& averaging(eslocal corner)
	{
		return _averaging[corner];
	}

	const std::vector<std::vector<eslocal> > boundary() const
	{
		return _boundaries;
	}

private:
	/** @brief Keeps mapping of nodes to mesh parts. */
	std::vector<std::vector<eslocal> > _boundaries;

	/** @brief Keeps information whether a point is the corner. */
	std::vector<bool> _corners;

	/** @brief Defines corners averaging nodes. */
	std::map<eslocal, std::vector<eslocal> > _averaging;
};

}

#endif /* BOUNDARIES_H_ */
