#ifndef SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_
#define SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_

#include "volumeelement.h"
#include "../plane/triangle3.h"
#include "../plane/triangle6.h"
#include "tetrahedron4.h"

#define Tetrahedron10NodesCount 10
#define Tetrahedron10EdgeCount 6
#define Tetrahedron10FacesCount 4
#define Tetrahedron10GPCount 15
#define Tetrahedron10CommonNodes 4
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

	virtual Point faceNormal(const Element *face) const;
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) const;
	virtual Element* face(size_t index) const { return _faces[index]; }
	virtual Element* edge(size_t index) const { return _edges[index]; }

	const std::vector<DenseMatrix>& dN() const { return Tetrahedron10::_dN; }
	const std::vector<DenseMatrix>& N() const { return Tetrahedron10::_N; }
	const std::vector<double>& weighFactor() const { return Tetrahedron10::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Tetrahedron10::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Tetrahedron10::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Tetrahedron10::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Tetrahedron10::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Tetrahedron10::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	void setFace(size_t index, Element* face) { _faces[index] = face; }
	void setEdge(size_t index, Element* edge) { _edges[index] = edge; }
	void setFace(Element* face);
	void setEdge(Element* edge);

	void fillFaces();
	void fillEdges();

private:
	eslocal _indices[Tetrahedron10NodesCount];

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

#endif /* SRC_MESH_ELEMENTS_VOLUME_TETRAHEDRON10_H_ */
