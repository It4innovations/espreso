
#ifndef TETRAHEDRON4_H_
#define TETRAHEDRON4_H_

#include "../plane/triangle3.h"
#include "../element.h"

#define Tetrahedron4NodesCount 4
#define Tetrahedron4EdgeCount 6
#define Tetrahedron4FacesCount 4
#define Tetrahedron4GPCount 4
#define Tetrahedron4CommonNodes 3
#define Tetrahedron4VTKCode 10

namespace espreso {

class Tetrahedron4: public Element
{
	friend class Tetrahedron10;

public:
	static bool match(const eslocal *indices, eslocal n);
	static size_t counter()
	{
		return _counter;
	}
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
	eslocal param(Params param) const { return _params[param]; };
	void setParam(Params param, eslocal value) { _params[param] = value; }
	size_t params() const { return PARAMS_SIZE; }

	size_t faces() const { return Tetrahedron4FacesCount; }
	size_t edges() const { return Tetrahedron4EdgeCount; }
	size_t nodes() const { return Tetrahedron4NodesCount; }
	size_t coarseNodes() const { return Tetrahedron4NodesCount; }
	size_t gaussePoints() const { return Tetrahedron4GPCount; }

	virtual Element* face(size_t index) const { return _faces[index]; }
	virtual Element* edge(size_t index) const { return _edges[index]; }

	const std::vector<DenseMatrix>& dN() const { return Tetrahedron4::_dN; }
	const std::vector<DenseMatrix>& N() const { return Tetrahedron4::_N; }
	const std::vector<double>& weighFactor() const { return Tetrahedron4::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Tetrahedron4::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Tetrahedron4::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Tetrahedron4::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Tetrahedron4::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Tetrahedron4::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	void setFace(size_t index, Element* face) { _faces[index] = face; }
	void setEdge(size_t index, Element* edge) { _edges[index] = edge; }

	void fillFaces();
	void fillEdges();

private:
	eslocal _indices[Tetrahedron4NodesCount];
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

#endif /* TETRAHEDRON4_H_ */
