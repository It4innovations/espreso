#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <cstring>

#include "1D/point.h"

#include "../structures/coordinates.h"
#include "esbasis.h"

namespace espreso {



class Element
{
public:
	enum IndicesType {
		GLOBAL,
		LOCAL
	};

	enum Params {
		MATERIAL,
		TYPE,
		CONSTANT,
		COORDINATES,
		BODY,
		NUMBER,
		PARAMS_SIZE
	};

	inline static bool match(const eslocal *indices, eslocal x, eslocal y)
	{
		return indices[x] == indices[y];
	}

	friend std::ostream& operator<<(std::ostream& os, const Element &e);
	friend std::ofstream& operator<<(std::ofstream& os, const Element &e);

	Element(const eslocal *params)
	{
		memcpy(_params, params, sizeof(eslocal) * PARAMS_SIZE);
	}

	Element(std::ifstream &is)
	{
		// TODO:
		for (size_t p = 0; p < PARAMS_SIZE; p++) {
			_params[p] = 0;
		}
	}

	virtual ~Element() {};

	const eslocal& node(size_t index) const
	{
		return indices()[index];
	}

	eslocal& node(size_t index)
	{
		return indices()[index];
	}

	void setParams(eslocal *params)
	{
		memcpy(_params, params, sizeof(eslocal) * PARAMS_SIZE);
	}

	const eslocal* getParams()
	{
		return _params;
	}

	void setParam(Params param, eslocal value)
	{
		_params[param] = value;
	}

	eslocal getParam(Params param)
	{
		return _params[param];
	}

	virtual Element* copy() const = 0;

	virtual const std::vector<DenseMatrix>& dN() const = 0;
	virtual const std::vector<DenseMatrix>&  N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	// Virtual methods
	virtual eslocal nCommon() const = 0;
	virtual eslocal vtkCode() const = 0;
	virtual size_t size() const = 0;
	virtual size_t coarseSize() const = 0;
	virtual size_t gpSize() const = 0;
	virtual size_t faces() const = 0;

	virtual std::vector<eslocal> getFace(size_t face) const = 0;
	virtual Element* getFullFace(size_t face) const = 0;
	virtual Element* getCoarseFace(size_t face) const = 0;


	virtual std::vector<eslocal> getNeighbours(size_t nodeIndex) const = 0;
	virtual const eslocal* indices() const = 0;

protected:
	virtual eslocal* indices() = 0;

	eslocal _params[PARAMS_SIZE];

};

}


#endif /* ELEMENT_H_ */
