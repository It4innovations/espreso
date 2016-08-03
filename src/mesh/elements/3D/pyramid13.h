
#ifndef PYRAMID13_H_
#define PYRAMID13_H_

#include "../element.h"
#include "pyramid5.h"
#include "../2D/triangle6.h"
#include "../2D/square8.h"

#define Pyramid13NodesCount 13
#define Pyramid13FacesCount 5
#define Pyramid13GPCount 8
#define Pyramid13VTKCode 27

namespace espreso {

class Pyramid13: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);
	static size_t counter()
	{
		return _counter;
	}
	static void setDOFs(
			const std::vector<Property> element,
			const std::vector<Property> face,
			const std::vector<Property> edge,
			const std::vector<Property> point,
			const std::vector<Property> midPoint)
	{
		_DOFElement = element;
		_DOFFace = face;
		_DOFEdge = edge;
		_DOFPoint = point;
		_DOFMidPoint = midPoint;
	}

	Pyramid13(const eslocal *indices, eslocal n, const eslocal *params);
	Pyramid13(std::ifstream &is);

	Element* copy() const
	{
		return new Pyramid13(*this);
	}

	eslocal vtkCode() const
	{
		return Pyramid13VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Pyramid13NodesCount;
	}

	size_t coarseSize() const
	{
		return Pyramid5NodesCount;
	}

	size_t gpSize() const
	{
		return Pyramid13GPCount;
	}

	size_t faces() const
	{
		return Pyramid13FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Pyramid13::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Pyramid13::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Pyramid13::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Pyramid13::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Pyramid13::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Pyramid13::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Pyramid13::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Pyramid13::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 3;
	}

	Element* getCoarseFace(size_t face) const
	{
		return Pyramid5::getF(_indices, _params, face);
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;
	Element* getFullFace(size_t face) const;

protected:

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Pyramid13NodesCount];

	static size_t _counter;

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;
};

}


#endif /* PYRAMID5_H_ */
