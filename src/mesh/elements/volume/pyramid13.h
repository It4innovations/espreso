
#ifndef SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_
#define SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_

#include "volumeelement.h"
#include "../plane/square8.h"
#include "../plane/triangle6.h"
#include "pyramid5.h"

#define Pyramid13NodesCount 13
#define Pyramid13EdgeCount 8
#define Pyramid13FacesCount 5
#define Pyramid13GPCount 14
#define Pyramid13CommonNodes 3
#define Pyramid13VTKCode 27

namespace espreso {

class Pyramid13: public VolumeElement
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

	Pyramid13(const eslocal *indices, eslocal n, const eslocal *params);
	Pyramid13(std::ifstream &is);
	Element* copy() const { return new Pyramid13(*this); }

	eslocal nCommon() const { return Pyramid13CommonNodes; }
	eslocal vtkCode() const { return Pyramid13VTKCode; }

	size_t faces() const { return Pyramid13FacesCount; }
	size_t edges() const { return Pyramid13EdgeCount; }
	size_t nodes() const { return Pyramid13NodesCount; }
	size_t coarseNodes() const { return Pyramid5NodesCount; }
	size_t gaussePoints() const { return Pyramid13GPCount; }

	virtual Point faceNormal(const Element *face) const;
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) const;

	const std::vector<DenseMatrix>& dN() const { return Pyramid13::_dN; }
	const std::vector<DenseMatrix>& N() const { return Pyramid13::_N; }
	const std::vector<double>& weighFactor() const { return Pyramid13::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Pyramid13::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Pyramid13::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Pyramid13::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Pyramid13::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Pyramid13::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	void fillFaces();
	void fillEdges();

private:
	eslocal _indices[Pyramid13NodesCount];

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


#endif /* SRC_MESH_ELEMENTS_VOLUME_PYRAMID13_H_ */
