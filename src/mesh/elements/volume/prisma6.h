
#ifndef SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_
#define SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_

#include "volumeelement.h"
#include "../line/line2.h"
#include "../plane/triangle3.h"
#include "../plane/square4.h"

#define Prisma6NodesCount 6
#define Prisma6EdgeCount 9
#define Prisma6FacesCount 5
#define Prisma6GPCount 9
#define Prisma6CommonNodes 3
#define Prisma6VTKCode 13

namespace espreso {

class Prisma6: public VolumeElement
{
	friend class Prisma15;

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

	Prisma6(const eslocal *indices, eslocal n, const eslocal *params);
	Prisma6(std::ifstream &is);
	Element* copy() const { return new Prisma6(*this); }

	eslocal nCommon() const { return Prisma6CommonNodes; }
	eslocal vtkCode() const { return Prisma6VTKCode; }

	size_t faces() const { return Prisma6FacesCount; }
	size_t edges() const { return Prisma6EdgeCount; }
	size_t nodes() const { return Prisma6NodesCount; }
	size_t coarseNodes() const { return Prisma6NodesCount; }
	size_t gaussePoints() const { return Prisma6GPCount; }

	const std::vector<eslocal>& faceNodes(size_t index) const { return Prisma6::_facesNodes[index]; }
	const std::vector<eslocal>& edgeNodes(size_t index) const { return Prisma6::_edgesNodes[index]; }

	const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return index < 3 ? Square4::_dN : Triangle3::_dN; }
	const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return index < 3 ? Square4::_N : Triangle3::_N; }
	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_dN; }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Line2::_N; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Prisma6::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Prisma6::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Prisma6::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Prisma6::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Prisma6::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Prisma6::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Prisma6::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Prisma6::_DOFMidPoint; }

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
	eslocal _indices[Prisma6NodesCount];

	static std::vector<Property> _DOFElement;
	static std::vector<Property> _DOFFace;
	static std::vector<Property> _DOFEdge;
	static std::vector<Property> _DOFPoint;
	static std::vector<Property> _DOFMidPoint;

	static std::vector<std::vector<eslocal> > _facesNodes;
	static std::vector<std::vector<eslocal> > _edgesNodes;
};

}


#endif /* SRC_MESH_ELEMENTS_VOLUME_PRISMA6_H_ */
