#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <cstring>

#include "1D/point.h"

#include "../structures/coordinates.h"
#include "../matrices/matrices.h"

namespace mesh {

class Element
{
public:
	enum IndicesType {
		GLOBAL,
		LOCAL
	};

	inline static bool match(eslocal *indices, eslocal x, eslocal y)
	{
		return indices[x] == indices[y];
	}

	friend std::ostream& operator<<(std::ostream& os, const Element &e);

	virtual ~Element() {};

	eslocal node(size_t index) const
	{
		return indices()[index];
	}

	void fillNodes(eslocal *nodes) const
	{
		memcpy(nodes, indices(), size() * sizeof(eslocal));
	}

	void setLocalIndices(std::vector<eslocal> &mapping)
	{
		for (size_t i = 0; i < size(); i++) {
			indices()[i] = mapping[node(i)];
		}
	}

	virtual Element* copy() const = 0;

	virtual const std::vector<DenseMatrix>& dN() const = 0;
	virtual const std::vector<DenseMatrix>&  N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	// Virtual methods
	virtual eslocal vtkCode() const = 0;
	virtual size_t size() const = 0;
	virtual size_t gpSize() const = 0;
	virtual size_t faces() const = 0;
	virtual std::vector<eslocal> getFace(size_t face) const = 0;
	virtual std::vector<eslocal> getNeighbours(size_t nodeIndex) const = 0;
	virtual const eslocal* indices() const = 0;

protected:
	virtual eslocal* indices() = 0;

};

}


#endif /* ELEMENT_H_ */
