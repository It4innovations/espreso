#ifndef SQUARE8_H_
#define SQUARE8_H_

#include "../element.h"
#include "../line/line2.h"
#include "square4.h"

#define Square8NodesCount 8
#define Square8FacesCount 4
#define Square8GPCount 9
#define Square8VTKCode 23

namespace espreso {

class Square8: public Element
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

	Square8(const eslocal *indices, const eslocal *params);

	Element* copy() const
	{
		return new Square8(*this);
	}

	eslocal vtkCode() const
	{
		return Square8VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Square8NodesCount;
	}

	size_t coarseSize() const
	{
		return Square4NodesCount;
	}

	size_t gpSize() const
	{
		return Square8GPCount;
	}

	size_t faces() const
	{
		return Square8FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Square8::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Square8::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Square8::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Square8::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Square8::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Square8::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Square8::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Square8::_DOFMidPoint;
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

	eslocal* indices()
	{
		return _indices;
	}

private:
	eslocal _indices[Square8NodesCount];

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


#endif /* SQUARE8_H_ */
