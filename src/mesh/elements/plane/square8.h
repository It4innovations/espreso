
#ifndef SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_
#define SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_

#include "planeelement.h"
#include "square4.h"
#include "../line/line3.h"

#define Square8NodesCount 8
#define Square8EdgeCount 4
#define Square8GPCount 9
#define Square8CommonNodes 2
#define Square8VTKCode 23

namespace espreso {

class Square8: public PlaneElement
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

	Square8(const eslocal *indices);
	Square8(const eslocal *indices, const eslocal *params);
	Square8(std::ifstream &is);
	Element* copy() const { return new Square8(*this); }

	std::vector<std::vector<eslocal> > triangularize() const;

	eslocal nCommon() const { return Square8CommonNodes; }
	eslocal vtkCode() const { return Square8VTKCode; }

	size_t edges() const { return Square8EdgeCount; }
	size_t nodes() const { return Square8NodesCount; }
	size_t coarseNodes() const { return Square4NodesCount; }
	size_t gaussePoints() const { return Square8GPCount; }

	const std::vector<eslocal>& edgeNodes(size_t index) const { return Square8::_edgesNodes[index]; }

	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_dN; }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_N; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Square8::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Square8::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Square8::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Square8::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Square8::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Square8::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Square8::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Square8::_DOFMidPoint; }

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	size_t fillEdges();

private:
	eslocal _indices[Square8NodesCount];

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;

	static std::vector<std::vector<eslocal> > _edgesNodes;
};

}


#endif /* SRC_MESH_ELEMENTS_PLANE_SQUARE8_H_ */
