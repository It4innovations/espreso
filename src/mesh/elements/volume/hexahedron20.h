
#ifndef HEXAHEDRON20_H_
#define HEXAHEDRON20_H_

#include "../element.h"
#include "../plane/square8.h"
#include "hexahedron8.h"

#define Hexahedron20NodesCount 20
#define Hexahedron20EdgeCount 12
#define Hexahedron20FacesCount 6
#define Hexahedron20GPCount 8
#define Hexahedron20CommonNodes 4
#define Hexahedron20VTKCode 25

namespace espreso {

class Hexahedron20: public Element
{

public:
	static bool match(const eslocal *indices, eslocal n);
	static size_t counter() { return _counter; }
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
	eslocal param(Params param) const { return _params[param]; }
	void setParam(Params param, eslocal value) { _params[param] = value; }
	size_t params() const { return PARAMS_SIZE; }

	size_t faces() const { return Hexahedron20FacesCount; }
	size_t edges() const { return Hexahedron20EdgeCount; }
	size_t nodes() const { return Hexahedron20NodesCount; }
	size_t coarseNodes() const { return Hexahedron8NodesCount; }
	size_t gaussePoints() const { return Hexahedron20GPCount; }

	virtual Point faceNormal(const Element *face);
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates);
	virtual Element* face(size_t index) const { return _faces[index]; }
	virtual Element* edge(size_t index) const { return _edges[index]; }

	const std::vector<DenseMatrix>& dN() const { return Hexahedron20::_dN; }
	const std::vector<DenseMatrix>& N() const { return Hexahedron20::_N; }
	const std::vector<double>& weighFactor() const { return Hexahedron20::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Hexahedron20::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Hexahedron20::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Hexahedron20::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Hexahedron20::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Hexahedron20::_DOFMidPoint; }

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
	eslocal _indices[Hexahedron20NodesCount];
	eslocal _params[PARAMS_SIZE];
	std::vector<Element*> _edges;
	std::vector<Element*> _faces;

	static size_t _counter;

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

#endif /* HEXAHEDRON20_H_ */
