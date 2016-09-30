
#ifndef SRC_MESH_ELEMENTS_LINE_LINE3_H_
#define SRC_MESH_ELEMENTS_LINE_LINE3_H_

#include "lineelement.h"
#include "line2.h"

#define Line3NodesCount 3
#define Line3GPCount 3
#define Line3CommonNodes 1
#define Line3VTKCode 4

namespace espreso {

class Line3: public LineElement
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

	Line3(const eslocal *indices);
	Line3(std::ifstream &is);
	Element* copy() const { return new Line3(*this); }

	eslocal nCommon() const { return Line3CommonNodes; }
	eslocal vtkCode() const { return Line3VTKCode; }

	size_t nodes() const { return Line3NodesCount; }
	size_t coarseNodes() const { return Line2NodesCount; }
	size_t gaussePoints() const { return Line3GPCount; }

	const std::vector<DenseMatrix>& dN() const { return Line3::_dN; }
	const std::vector<DenseMatrix>& N() const { return Line3::_N; }
	const std::vector<double>& weighFactor() const { return Line3::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Line3::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Line3::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Line3::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Line3::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Line3::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

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


#endif /* SRC_MESH_ELEMENTS_LINE_LINE3_H_ */
