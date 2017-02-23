#ifndef SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_
#define SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_

#include "volumeelement.h"
#include "../line/line3.h"
#include "../plane/triangle3.h"
#include "../plane/triangle6.h"
#include "tetrahedron4.h"

#define Tetrahedron10NodesCount 10
#define Tetrahedron10EdgeCount 6
#define Tetrahedron10FacesCount 4
// WARNINKG: use only 15 gausse points
#define Tetrahedron10GPCount 15
#define Tetrahedron10CommonNodes 3
#define Tetrahedron10VTKCode 24

namespace espreso {

class Tetrahedron10: public VolumeElement
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

	Tetrahedron10(const eslocal *indices, eslocal n, const eslocal *params);
	Tetrahedron10(std::ifstream &is);
	Element* copy() const { return new Tetrahedron10(*this); }

	eslocal nCommon() const { return Tetrahedron10CommonNodes; }
	eslocal vtkCode() const { return Tetrahedron10VTKCode; }

	size_t faces() const { return Tetrahedron10FacesCount; }
	size_t edges() const { return Tetrahedron10EdgeCount; }
	size_t nodes() const { return Tetrahedron10NodesCount; }
	size_t coarseNodes() const { return Tetrahedron4NodesCount; }
	size_t gaussePoints() const { return Tetrahedron10GPCount; }

	const std::vector<eslocal>& faceNodes(size_t index) const { return Tetrahedron10::_facesNodes[index]; }
	const std::vector<eslocal>& edgeNodes(size_t index) const { return Tetrahedron10::_edgesNodes[index]; }

	const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Triangle6::_dN; }
	const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Triangle6::_N; }
	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_dN; }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line3::_N; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Tetrahedron10::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Tetrahedron10::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Tetrahedron10::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Tetrahedron10::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Tetrahedron10::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Tetrahedron10::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Tetrahedron10::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Tetrahedron10::_DOFMidPoint; }

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
	eslocal _indices[Tetrahedron10NodesCount];

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;

	static std::vector<std::vector<eslocal> > _facesNodes;
	static std::vector<std::vector<eslocal> > _edgesNodes;
};

}

#endif /* SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_ */
