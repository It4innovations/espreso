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

	inline static bool match(esint *indices, esint x, esint y)
	{
		return indices[x] == indices[y];
	}

	friend std::ostream& operator<<(std::ostream& os, const Element &e);

	virtual ~Element() {};

	esint node(size_t index) const
	{
		return indices()[index];
	}

	void fillNodes(esint *nodes) const
	{
		memcpy(nodes, indices(), size() * sizeof(esint));
	}

	void setLocalIndices(std::vector<esint> &mapping)
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
	virtual esint vtkCode() const = 0;
	virtual size_t size() const = 0;
	virtual size_t gpSize() const = 0;
	virtual size_t faces() const = 0;
	virtual std::vector<esint> getFace(size_t face) const = 0;
	virtual std::vector<esint> getNeighbours(size_t nodeIndex) const = 0;
	virtual const esint* indices() const = 0;

protected:
	virtual esint* indices() = 0;

};

}


#endif /* ELEMENT_H_ */
