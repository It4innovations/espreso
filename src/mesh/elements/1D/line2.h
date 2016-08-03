#ifndef LINE2_H_
#define LINE2_H_


#include "../element.h"

#define Line2NodesCount 2
#define Line2FacesCount 0
#define Line2GPCount 2
#define Line2VTKCode 3

namespace espreso {

class Line2: public Element
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

	Line2(const eslocal *indices, const eslocal *params);

	Element* copy() const
	{
		return new Line2(*this);
	}

	eslocal vtkCode() const
	{
		return Line2VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Line2NodesCount;
	}

	size_t coarseSize() const
	{
		return Line2NodesCount;
	}

	size_t gpSize() const
	{
		return Line2GPCount;
	}

	size_t faces() const
	{
		return Line2FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Line2::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Line2::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Line2::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Line2::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Line2::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Line2::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Line2::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Line2::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 1;
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
	eslocal _indices[Line2NodesCount];

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


#endif /* LINE2_H_ */
