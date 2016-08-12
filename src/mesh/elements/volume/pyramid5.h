
#ifndef PYRAMID5_H_
#define PYRAMID5_H_

#include "../element.h"
#include "../plane/triangle3.h"
#include "../plane/square4.h"

#define Pyramid5NodesCount 5
#define Pyramid5EdgeCount 8
#define Pyramid5FacesCount 5
#define Pyramid5GPCount 8
#define Pyramid5CommonNodes 3
#define Pyramid5VTKCode 14

namespace espreso {

class Pyramid5: public Element
{
	friend class Pyramid13;

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

	Pyramid5(const eslocal *indices, eslocal n, const eslocal *params);
	Pyramid5(std::ifstream &is);
	Element* copy() const { return new Pyramid5(*this); }

	eslocal nCommon() const { return Pyramid5CommonNodes; }
	eslocal vtkCode() const { return Pyramid5VTKCode; }
	eslocal param(Params param) const { return _params[param]; };
	void param(Params param, eslocal value) { _params[param] = value; }

	size_t faces() const { return Pyramid5FacesCount; }
	size_t edges() const { return Pyramid5EdgeCount; }
	size_t nodes() const { return Pyramid5NodesCount; }
	size_t coarseNodes() const { return Pyramid5NodesCount; }
	size_t gaussePoints() const { return Pyramid5GPCount; }

	virtual Element* face(size_t index) const { return _faces[index]; }
	virtual Element* edge(size_t index) const { return _edges[index]; }

	const std::vector<DenseMatrix>& dN() const { return Pyramid5::_dN; }
	const std::vector<DenseMatrix>& N() const { return Pyramid5::_N; }
	const std::vector<double>& weighFactor() const { return Pyramid5::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Pyramid5::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Pyramid5::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Pyramid5::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Pyramid5::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Pyramid5::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const;
	eslocal* indices() { return _indices; }
	const eslocal* indices() const { return _indices; }

	void setFace(size_t index, Element* face) { _faces[index] = face; }
	void setEdge(size_t index, Element* edge) { _edges[index] = edge; }

	void fillFaces();
	void fillEdges();

private:
	eslocal _indices[Pyramid5NodesCount];
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


#endif /* PYRAMID5_H_ */
