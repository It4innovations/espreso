#ifndef TRIANGLE3_H_
#define TRIANGLE3_H_

#include "../element.h"
#include "../line/line2.h"

#define Triangle3NodesCount 3
#define Triangle3FacesCount 3
#define Triangle3GPCount 1
#define Triangle3VTKCode 5

namespace espreso {

class Triangle3: public Element
{

public:
	static bool match(eslocal *indices, eslocal n);
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

	Triangle3(const eslocal *indices, const eslocal *params);

	Element* copy() const
	{
		return new Triangle3(*this);
	}

	eslocal vtkCode() const
	{
		return Triangle3VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Triangle3NodesCount;
	}

	size_t coarseSize() const
	{
		return Triangle3NodesCount;
	}

	size_t gpSize() const
	{
		return Triangle3GPCount;
	}

	size_t faces() const
	{
		return Triangle3FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Triangle3::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Triangle3::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Triangle3::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Triangle3::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Triangle3::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Triangle3::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Triangle3::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Triangle3::_DOFMidPoint;
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
	eslocal _indices[Triangle3NodesCount];

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



#endif /* TRIANGLE3_H_ */
