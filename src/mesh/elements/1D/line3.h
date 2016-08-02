#ifndef LINE3_H_
#define LINE3_H_

#include "../element.h"
#include "line2.h"

#define Line3NodesCount 3
#define Line3FacesCount 0
#define Line3GPCount 2
#define Line3VTKCode 4

namespace espreso {

class Line3: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);

	Line3(const eslocal *indices, const eslocal *params);

	Element* copy() const
	{
		return new Line3(*this);
	}

	eslocal vtkCode() const
	{
		return Line3VTKCode;
	}

	const eslocal* indices() const
	{
		return _indices;
	}

	size_t size() const
	{
		return Line3NodesCount;
	}

	size_t coarseSize() const
	{
		return Line2NodesCount;
	}

	size_t gpSize() const
	{
		return Line3GPCount;
	}

	size_t faces() const
	{
		return Line3FacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Line3::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Line3::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Line3::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Line3::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Line3::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Line3::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Line3::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Line3::_DOFMidPoint;
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
	eslocal _indices[Line3NodesCount];

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


#endif /* LINE3_H_ */
