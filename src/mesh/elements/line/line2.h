
#ifndef SRC_MESH_ELEMENTS_LINE_LINE2_H_
#define SRC_MESH_ELEMENTS_LINE_LINE2_H_

#include "lineelement.h"

#define Line2NodesCount 2
#define Line2GPCount 2
#define Line2CommonNodes 1
#define Line2VTKCode 3

namespace espreso {

class Line2: public LineElement
{

public:
	static bool match(const eslocal *indices, eslocal n);
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

	Line2(const eslocal *indices);
	Line2(std::ifstream &is);
	Element* copy() const { return new Line2(*this); }

	eslocal nCommon() const { return Line2CommonNodes; }
	eslocal vtkCode() const { return Line2VTKCode; }

	size_t nodes() const { return Line2NodesCount; }
	size_t coarseNodes() const { return Line2NodesCount; }
	size_t gaussePoints() const { return Line2GPCount; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Line2::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Line2::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Line2::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Line2::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Line2::_DOFMidPoint; }


protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	const eslocal* indices() const { return _indices; }
	eslocal* indices() { return _indices; }

private:
	eslocal _indices[Line2NodesCount];

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


#endif /* SRC_MESH_ELEMENTS_LINE_LINE2_H_ */
