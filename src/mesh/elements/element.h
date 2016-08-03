#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <cstring>

#include "../structures/coordinates.h"
#include "esbasis.h"

namespace espreso {

class Mesh;

class Element
{
	friend class Mesh;

public:
	enum IndicesType {
		GLOBAL,
		LOCAL
	};

	enum Params {
		MATERIAL,
		CONSTANT,
		COORDINATES,
		BODY,
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

	inline bool operator==(const Element& other)
	{
		if (size() != other.size()) {
			return false;
		}
		return std::is_permutation(indices(), indices() + size(), other.indices());
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

	const eslocal* getParams() const
	{
		return _params;
	}

	void setParam(Params param, eslocal value)
	{
		_params[param] = value;
	}

	eslocal getParam(Params param) const
	{
		return _params[param];
	}

	virtual Element* copy() const = 0;

	virtual const std::vector<DenseMatrix>& dN() const = 0;
	virtual const std::vector<DenseMatrix>&  N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	virtual const std::vector<Property>& DOFElement() const = 0;
	virtual const std::vector<Property>& DOFFace() const = 0;
	virtual const std::vector<Property>& DOFEdge() const = 0;
	virtual const std::vector<Property>& DOFPoint() const = 0;
	virtual const std::vector<Property>& DOFMidPoint() const = 0;

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


	Settings& settings()
	{
		return _settings;
	}

	const Settings& settings() const
	{
		return _settings;
	}

	std::vector<Evaluator*>& settings(Property property)
	{
		return _settings[property];
	}

	const std::vector<Evaluator*>& settings(Property property) const
	{
		return _settings[property];
	}

	/// Returns sub-domains where is the element
	const std::vector<eslocal>& domains() const
	{
		return _domains;
	}

	/// Returns clusters where is the element
	const std::vector<eslocal>& clusters() const
	{
		return _clusters;
	}

	std::vector<eslocal>& clusters()
	{
		return _clusters;
	}

protected:
	virtual eslocal* indices() = 0;

	eslocal _params[PARAMS_SIZE];

	Settings _settings;
	std::vector<eslocal> _domains;
	std::vector<eslocal> _clusters;

};

}


#endif /* ELEMENT_H_ */
