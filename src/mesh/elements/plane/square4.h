#ifndef SQUARE4_H_
#define SQUARE4_H_

#include "../element.h"
#include "../line/line2.h"

#define Square4NodesCount 4
#define Square4EdgeCount 4
#define Square4FacesCount 4
#define Square4GPCount 4
#define Square4CommonNodes 2
#define Square4VTKCode 9

namespace espreso {

class Square4: public Element
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

	Square4(const eslocal *indices);
	Square4(const eslocal *indices, const eslocal *params);
	Element* copy() const { return new Square4(*this); }

	eslocal nCommon() const { return Square4CommonNodes; }
	eslocal vtkCode() const { return Square4VTKCode; }
	eslocal param(Params param) const { return _params[param]; }
	void param(Params param, eslocal value) { _params[param] = value; }

	size_t faces() const { return Square4FacesCount; }
	size_t edges() const { return Square4EdgeCount; }
	size_t nodes() const { return Square4NodesCount; }
	size_t coarseNodes() const { return Square4NodesCount; }
	size_t gaussePoints() const { return Square4GPCount; }

	virtual Element* face(size_t index) const { ESINFO(GLOBAL_ERROR) << "Square4 has no face"; return NULL; }
	virtual Element* edge(size_t index) const { return _edges[index]; };

	const std::vector<DenseMatrix>& dN() const { return Square4::_dN; }
	const std::vector<DenseMatrix>& N() const { return Square4::_N; }
	const std::vector<double>& weighFactor() const { return Square4::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Square4::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Square4::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Square4::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Square4::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Square4::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	void face(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Square4 has no face"; }
	void edge(size_t index, Element* edge) { _edges[index] = edge; }

	void fillFaces() {};
	void fillEdges();

private:
	eslocal _indices[Square4NodesCount];
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


#endif /* SQUARE4_H_ */
