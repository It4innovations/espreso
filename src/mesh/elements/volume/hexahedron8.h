#ifndef HEXAHEDRON8_H_
#define HEXAHEDRON8_H_

#include "../element.h"
#include "../plane/square4.h"
#include "../line/line2.h"

#define Hexahedron8NodesCount 8
#define Hexahedron8FacesCount 6
#define Hexahedron8GPCount 8
#define Hexahedron8VTKCode 12

namespace espreso {

class Hexahedron8: public Element
{
	friend class Hexahedron20;

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

	Hexahedron8(const eslocal *indices, eslocal n, const eslocal *params);
	Hexahedron8(std::ifstream &is);

	Element* copy() const
	{
		return new Hexahedron8(*this);
	}

	eslocal vtkCode() const
	{
		return Hexahedron8VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Hexahedron8NodesCount;
	}

	size_t coarseSize() const
	{
		return Hexahedron8NodesCount;
	}

	size_t gpSize() const
	{
		return Hexahedron8GPCount;
	}

	size_t faces() const
	{
		return Hexahedron8FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Hexahedron8::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Hexahedron8::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Hexahedron8::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Hexahedron8::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Hexahedron8::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Hexahedron8::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Hexahedron8::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Hexahedron8::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 3;
	}

	Element* getFullFace(size_t face) const
	{
		return getF(_indices, _params, face);
	}

	Element* getCoarseFace(size_t face) const
	{
		return getFullFace(face);
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;


protected:
	static Element* getF(const eslocal *indices, const eslocal *params, size_t face);

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Hexahedron8NodesCount];

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

#endif /* HEXAHEDRON8_H_ */
