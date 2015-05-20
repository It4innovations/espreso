#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <cstring>

#include "1D/point.h"

#include "../structures/coordinates.h"
#include "../matrices/matrices.h"

class Element
{
public:
	enum IndicesType {
		GLOBAL,
		LOCAL
	};

	inline static bool match(idx_t *indices, idx_t x, idx_t y)
	{
		return indices[x] == indices[y];
	}

	friend std::ostream& operator<<(std::ostream& os, const Element &e);

	virtual ~Element() {};

	idx_t node(size_t index) const
	{
		return indices()[index];
	}

	void fillNodes(idx_t *nodes) const
	{
		for (size_t i = 0; i < size(); i++) {
			nodes[i] = node(i);
		}
	}
	void setLocalIndices(std::vector<idx_t> &mapping)
	{
		for (size_t i = 0; i < size(); i++) {
			indices()[i] = mapping[node(i)];
		}
	}

	virtual Element* copy() const = 0;

	virtual const std::vector<std::vector<double> >& dN() const = 0;
	virtual const std::vector<std::vector<double> >&  N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	// Virtual methods
	virtual int vtkCode() const = 0;
	virtual size_t size() const = 0;
	virtual size_t gpSize() const = 0;
	virtual size_t faces() const = 0;
	virtual std::vector<idx_t> getFace(size_t face) const = 0;
	virtual std::vector<idx_t> getNeighbours(size_t nodeIndex) const = 0;
	virtual const idx_t* indices() const = 0;

protected:
	virtual idx_t* indices() = 0;

};




#endif /* ELEMENT_H_ */
