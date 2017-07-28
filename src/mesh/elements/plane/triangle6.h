
#ifndef SRC_MESH_ELEMENTS_PLANE_TRIANGLE6_H_
#define SRC_MESH_ELEMENTS_PLANE_TRIANGLE6_H_

#include "planeelement.h"
#include "triangle3.h"
#include "../line/line3.h"

#define Triangle6NodesCount 6
#define Triangle6EdgeCount 3
#define Triangle6GPCount 6
#define Triangle6CommonNodes 2
#define Triangle6VTKCode 22

namespace espreso {

class Triangle6: public PlaneElement
{

public:
	static bool match(const eslocal *indices, const eslocal n);
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

	Triangle6(const eslocal *indices);
	Triangle6(const eslocal *indices, const eslocal *params);
	Triangle6(std::ifstream &is);
	Element* copy() const { return new Triangle6(*this); }

	std::vector<std::vector<eslocal> > triangularize() const;

	eslocal nCommon() const { return Triangle6CommonNodes; }
	eslocal vtkCode() const { return Triangle6VTKCode; }

	size_t edges() const { return Triangle6EdgeCount; }
	size_t nodes() const { return Triangle6NodesCount; }
	size_t coarseNodes() const { return Triangle3NodesCount; }
	size_t gaussePoints() const { return Triangle6GPCount; }

	const std::vector<eslocal>& edgeNodes(size_t index) const { return Triangle6::_edgesNodes[index]; }

	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_dN; }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_N; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Triangle6::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Triangle6::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Triangle6::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Triangle6::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Triangle6::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Triangle6::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Triangle6::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Triangle6::_DOFMidPoint; }

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	size_t fillEdges();

private:
	eslocal _indices[Triangle6NodesCount];

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;

	static std::vector<std::vector<eslocal> > _edgesNodes;
};

}



#endif /* SRC_MESH_ELEMENTS_PLANE_TRIANGLE6_H_ */
