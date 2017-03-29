
#ifndef SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_
#define SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_

#include "volumeelement.h"
#include "../plane/square4.h"
#include "../line/line2.h"

#define Hexahedron8NodesCount 8
#define Hexahedron8EdgeCount 12
#define Hexahedron8FacesCount 6
#define Hexahedron8GPCount 8
#define Hexahedron8CommonNodes 3
#define Hexahedron8VTKCode 12

namespace espreso {

class Hexahedron8: public VolumeElement
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

	Hexahedron8(const eslocal *indices, eslocal n, const eslocal *params);
	Hexahedron8(std::ifstream &is);
	Element* copy() const { return new Hexahedron8(*this); }

	eslocal nCommon() const { return Hexahedron8CommonNodes; }
	eslocal vtkCode() const { return Hexahedron8VTKCode; }

	size_t faces() const { return Hexahedron8FacesCount; }
	size_t edges() const { return Hexahedron8EdgeCount; }
	size_t nodes() const { return Hexahedron8NodesCount; }
	size_t coarseNodes() const { return Hexahedron8NodesCount; }
	size_t gaussePoints() const { return Hexahedron8GPCount; }

	Element* addFace(const std::vector<eslocal> &nodes);

	const std::vector<eslocal>& faceNodes(size_t index) const { return Hexahedron8::_facesNodes[index]; }
	const std::vector<eslocal>& edgeNodes(size_t index) const { return Hexahedron8::_edgesNodes[index]; }

	const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Square4::_dN; }
	const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Square4::_N; }
	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_dN; }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_N; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Hexahedron8::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Hexahedron8::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Hexahedron8::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Hexahedron8::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Hexahedron8::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Hexahedron8::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Hexahedron8::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Hexahedron8::_DOFMidPoint; }

	static std::vector<DenseMatrix> _dN;
	static std::vector<DenseMatrix> _N;
	static std::vector<double> _weighFactor;

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	size_t fillFaces();
	size_t fillEdges();

private:
	eslocal _indices[Hexahedron8NodesCount];

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;

	static std::vector<std::vector<eslocal> > _facesNodes;
	static std::vector<std::vector<eslocal> > _edgesNodes;
};

}

#endif /* SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON8_H_ */
