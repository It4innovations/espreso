
#ifndef SRC_MESH_ELEMENTS_POINT_NODE_H_
#define SRC_MESH_ELEMENTS_POINT_NODE_H_

#include "../element.h"

#define NodeNodesCount 1
#define NodeEdgeCount 0
#define NodeFacesCount 0
#define NodeGPCount 0
#define NodeCommonNodes 1
#define NodeVTKCode 2

namespace espreso {

class Node: public Element
{

public:
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

	Node(eslocal index): _index(index) {};
	Element* copy() const { return new Node(*this); }

	eslocal nCommon() const { return NodeCommonNodes; }
	eslocal vtkCode() const { return NodeVTKCode; }
	eslocal param(Params param) const { ESINFO(GLOBAL_ERROR) << "Node has no params"; return 0; }
	void param(Params param, eslocal value) { ESINFO(GLOBAL_ERROR) << "Node has no params"; }

	size_t faces() const { return NodeFacesCount; }
	size_t edges() const { return NodeEdgeCount; }
	size_t nodes() const { return NodeNodesCount; }
	size_t coarseNodes() const { return NodeNodesCount; }
	size_t gaussePoints() const { return NodeGPCount; }

	virtual Element* face(size_t index) const { ESINFO(GLOBAL_ERROR) << "Node has no face"; return NULL; }
	virtual Element* edge(size_t index) const { ESINFO(GLOBAL_ERROR) << "Node has no edge"; return NULL; }

	const std::vector<DenseMatrix>& dN() const { return Node::_dN; }
	const std::vector<DenseMatrix>& N() const { return Node::_N; }
	const std::vector<double>& weighFactor() const { return Node::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Node::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Node::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Node::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Node::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Node::_DOFMidPoint; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const { ESINFO(GLOBAL_ERROR) << "Node has no neighbours"; return {}; }
	eslocal* indices() { return &_index; }
	const eslocal* indices() const { return &_index; }

	void face(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Node has no face"; }
	void edge(size_t index, Element* edge) { ESINFO(GLOBAL_ERROR) << "Node has no edge"; }

	void fillFaces() {};
	void fillEdges() {};

private:
	eslocal _index;

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



#endif /* SRC_MESH_ELEMENTS_POINT_NODE_H_ */
