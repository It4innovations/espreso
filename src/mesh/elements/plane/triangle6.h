#ifndef TRIANGLE6_H_
#define TRIANGLE6_H_

#include "../element.h"
#include "triangle3.h"

#define Triangle6NodesCount 6
#define Triangle6EdgeCount 3
#define Triangle6FacesCount 0
#define Triangle6GPCount 6
#define Triangle6CommonNodes 3
#define Triangle6VTKCode 22

namespace espreso {

class Triangle6: public Element
{

public:
	static bool match(const eslocal *indices, const eslocal n);
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

	Triangle6(const eslocal *indices);
	Triangle6(const eslocal *indices, const eslocal *params);
	Element* copy() const { return new Triangle6(*this); }

	eslocal nCommon() const { return Triangle6CommonNodes; }
	eslocal vtkCode() const { return Triangle6VTKCode; }
	eslocal param(Params param) const { return _params[param]; };
	void setParam(Params param, eslocal value) { _params[param] = value; }
	size_t params() const { return _params.size(); }

	size_t faces() const { return Triangle6FacesCount; }
	size_t edges() const { return Triangle6EdgeCount; }
	size_t nodes() const { return Triangle6NodesCount; }
	size_t coarseNodes() const { return Triangle3NodesCount; }
	size_t gaussePoints() const { return Triangle6GPCount; }

	virtual Point faceNormal(const Element *face) { ESINFO(GLOBAL_ERROR) << "Triangle6 has no face"; return Point(); }
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates);
	virtual Element* face(size_t index) const { ESINFO(GLOBAL_ERROR) << "Triangle6 has no face"; return NULL; }
	virtual Element* edge(size_t index) const { return _edges[index]; };

	const std::vector<DenseMatrix>& dN() const { return Triangle6::_dN; }
	const std::vector<DenseMatrix>& N() const { return Triangle6::_N; }
	const std::vector<double>& weighFactor() const { return Triangle6::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Triangle6::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Triangle6::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Triangle6::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Triangle6::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Triangle6::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	void setFace(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Triangle6 has no face"; }
	void setEdge(size_t index, Element* edge) { _edges[index] = edge; }
	void setFace(Element* face) { ESINFO(GLOBAL_ERROR) << "Triangle6 has no face"; }
	void setEdge(Element* edge);

	void fillFaces() {};
	void fillEdges();

private:
	eslocal _indices[Triangle6NodesCount];
	std::vector<eslocal> _params;
	std::vector<Element*> _edges;

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



#endif /* TRIANGLE6_H_ */
