#ifndef LINE2_H_
#define LINE2_H_

#include "../element.h"

#define Line2NodesCount 2
#define Line2EdgeCount 0
#define Line2FacesCount 0
#define Line2GPCount 2
#define Line2CommonNodes 1
#define Line2VTKCode 3

namespace espreso {

class Line2: public Element
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

	Line2(const eslocal *indices);
	Line2(std::ifstream &is);
	Element* copy() const { return new Line2(*this); }

	eslocal nCommon() const { return Line2CommonNodes; }
	eslocal vtkCode() const { return Line2VTKCode; }
	eslocal param(Params param) const { ESINFO(GLOBAL_ERROR) << "Line2 has no params"; return 0; }
	void setParam(Params param, eslocal value) { ESINFO(GLOBAL_ERROR) << "Line2 has no params"; }
	size_t params() const { return 0; }

	size_t faces() const { return Line2FacesCount; }
	size_t edges() const { return Line2EdgeCount; }
	size_t nodes() const { return Line2NodesCount; }
	size_t coarseNodes() const { return Line2NodesCount; }
	size_t gaussePoints() const { return Line2GPCount; }

	virtual Point faceNormal(const Element *face) { ESINFO(GLOBAL_ERROR) << "Line2 has no face"; return Point(); }
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) { ESINFO(GLOBAL_ERROR) << "Line2 has no edge"; return Point(); }
	virtual Element* face(size_t index) const { ESINFO(GLOBAL_ERROR) << "Line2 has no face"; return NULL; }
	virtual Element* edge(size_t index) const { ESINFO(GLOBAL_ERROR) << "Line2 has no edge"; return NULL; }

	const std::vector<DenseMatrix>& dN() const { return Line2::_dN; }
	const std::vector<DenseMatrix>& N() const { return Line2::_N; }
	const std::vector<double>& weighFactor() const { return Line2::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Line2::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Line2::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Line2::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Line2::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Line2::_DOFMidPoint; }


protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	const eslocal* indices() const { return _indices; }
	eslocal* indices() { return _indices; }

	void setFace(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Line2 has no face"; }
	void setEdge(size_t index, Element* edge) { ESINFO(GLOBAL_ERROR) << "Line2 has no edge"; }
	void setFace(Element* face) { ESINFO(GLOBAL_ERROR) << "Line2 has no face"; }
	void setEdge(Element* edge) { ESINFO(GLOBAL_ERROR) << "Line2 has no edge"; }

	void fillFaces() {};
	void fillEdges() {};

private:
	eslocal _indices[Line2NodesCount];

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


#endif /* LINE2_H_ */
