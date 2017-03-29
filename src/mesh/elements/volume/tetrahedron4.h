
#ifndef SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_
#define SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_

#include "../line/line2.h"
#include "../plane/triangle3.h"
#include "volumeelement.h"

#define Tetrahedron4NodesCount 4
#define Tetrahedron4EdgeCount 6
#define Tetrahedron4FacesCount 4
#define Tetrahedron4GPCount 4
#define Tetrahedron4CommonNodes 3
#define Tetrahedron4VTKCode 10

namespace espreso {

class Tetrahedron4: public VolumeElement
{
	friend class Tetrahedron10;

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

	Tetrahedron4(const eslocal *indices, eslocal n, const eslocal *params);
	Tetrahedron4(std::ifstream &is);
	Element* copy() const { return new Tetrahedron4(*this); }

	eslocal nCommon() const { return Tetrahedron4CommonNodes; }
	eslocal vtkCode() const { return Tetrahedron4VTKCode; }

	size_t faces() const { return Tetrahedron4FacesCount; }
	size_t edges() const { return Tetrahedron4EdgeCount; }
	size_t nodes() const { return Tetrahedron4NodesCount; }
	size_t coarseNodes() const { return Tetrahedron4NodesCount; }
	size_t gaussePoints() const { return Tetrahedron4GPCount; }

	Element* addFace(const std::vector<eslocal> &nodes);

	const std::vector<eslocal>& faceNodes(size_t index) const { return Tetrahedron4::_facesNodes[index]; }
	const std::vector<eslocal>& edgeNodes(size_t index) const { return Tetrahedron4::_edgesNodes[index]; }

	const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Triangle3::_dN; }
	const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Triangle3::_N; }
	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_dN; }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_N; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Tetrahedron4::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Tetrahedron4::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Tetrahedron4::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Tetrahedron4::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Tetrahedron4::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Tetrahedron4::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Tetrahedron4::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Tetrahedron4::_DOFMidPoint; }

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
	eslocal _indices[Tetrahedron4NodesCount];

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;

	static std::vector<std::vector<eslocal> > _facesNodes;
	static std::vector<std::vector<eslocal> > _edgesNodes;
};


}

#endif /* SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON4_H_ */
