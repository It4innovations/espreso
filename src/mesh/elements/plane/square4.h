#ifndef SQUARE4_H_
#define SQUARE4_H_

#include "../element.h"
#include "../line/line2.h"

#define Square4NodesCount 4
#define Square4FacesCount 4
#define Square4GPCount 4
#define Square4VTKCode 9

namespace espreso {

class Square4: public Element
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

	Square4(const eslocal *indices, const eslocal *params);

	Element* copy() const
	{
		return new Square4(*this);
	}

	eslocal vtkCode() const
	{
		return Square4VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Square4NodesCount;
	}

	size_t coarseSize() const
	{
		return Square4NodesCount;
	}

	size_t gpSize() const
	{
		return Square4GPCount;
	}

	size_t faces() const
	{
		return Square4FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Square4::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Square4::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Square4::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Square4::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Square4::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Square4::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Square4::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Square4::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 2;
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	std::vector<eslocal> getFace(size_t face) const;
	Element* getFullFace(size_t face) const;
	Element* getCoarseFace(size_t face) const;

protected:
	static Element* getF(const eslocal *indices, const eslocal *params, size_t face);

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Square4NodesCount];

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


#endif /* SQUARE4_H_ */
