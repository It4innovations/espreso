
#ifndef SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_
#define SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_

#include "volumeelement.h"
#include "../plane/square8.h"
#include "hexahedron8.h"

#define Hexahedron20NodesCount 20
#define Hexahedron20EdgeCount 12
#define Hexahedron20FacesCount 6
#define Hexahedron20GPCount 14
#define Hexahedron20CommonNodes 4
#define Hexahedron20VTKCode 25

namespace espreso {

class Hexahedron20: public VolumeElement
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

	Hexahedron20(const eslocal *indices, eslocal n, const eslocal *params);
	Hexahedron20(std::ifstream &is);
	Element* copy() const { return new Hexahedron20(*this); }

	eslocal nCommon() const { return Hexahedron20CommonNodes; }
	eslocal vtkCode() const { return Hexahedron20VTKCode; }

	size_t faces() const { return Hexahedron20FacesCount; }
	size_t edges() const { return Hexahedron20EdgeCount; }
	size_t nodes() const { return Hexahedron20NodesCount; }
	size_t coarseNodes() const { return Hexahedron8NodesCount; }
	size_t gaussePoints() const { return Hexahedron20GPCount; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Hexahedron20::_dN; }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Hexahedron20::_N; }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { return Hexahedron20::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Hexahedron20::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Hexahedron20::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Hexahedron20::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Hexahedron20::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Hexahedron20::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	size_t fillFaces();
	size_t fillEdges();

private:
	eslocal _indices[Hexahedron20NodesCount];

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

#endif /* SRC_MESH_ELEMENTS_VOLUME_HEXAHEDRON20_H_ */
