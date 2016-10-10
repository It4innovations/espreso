
#ifndef SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_
#define SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_

#include "volumeelement.h"
#include "../plane/square8.h"
#include "../plane/triangle6.h"
#include "prisma6.h"

#define Prisma15NodesCount 15
#define Prisma15EdgeCount 9
#define Prisma15FacesCount 5
#define Prisma15GPCount 9
#define Prisma15CommonNodes 3
#define Prisma15VTKCode 26

namespace espreso {

class Prisma15: public VolumeElement
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

	Prisma15(const eslocal *indices, eslocal n, const eslocal *params);
	Prisma15(std::ifstream &is);
	Element* copy() const { return new Prisma15(*this); }

	eslocal nCommon() const { return Prisma15CommonNodes; }
	eslocal vtkCode() const { return Prisma15VTKCode; }

	size_t faces() const { return Prisma15FacesCount; }
	size_t edges() const { return Prisma15EdgeCount; }
	size_t nodes() const { return Prisma15NodesCount; }
	size_t coarseNodes() const { return Prisma6NodesCount; }
	size_t gaussePoints() const { return Prisma15GPCount; }

	virtual Point faceNormal(const Element *face) const;
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) const;
	virtual Element* face(size_t index) const { return _faces[index]; }
	virtual Element* edge(size_t index) const { return _edges[index]; }

	const std::vector<DenseMatrix>& dN() const { return Prisma15::_dN; }
	const std::vector<DenseMatrix>& N() const { return Prisma15::_N; }
	const std::vector<double>& weighFactor() const { return Prisma15::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Prisma15::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Prisma15::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Prisma15::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Prisma15::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Prisma15::_DOFMidPoint; }

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
	eslocal _indices[Prisma15NodesCount];

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


#endif /* SRC_MESH_ELEMENTS_VOLUME_PRISMA15_H_ */
