
#ifndef SRC_MESH_ELEMENTS_POINT_NODE_H_
#define SRC_MESH_ELEMENTS_POINT_NODE_H_

#include "pointelement.h"

#define NodeGPCount 0
#define NodeVTKCode 2

namespace espreso {

class Node: public PointElement
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

	Node(eslocal index): PointElement(index) {};
	Element* copy() const { return new Node(*this); }

	eslocal vtkCode() const { return NodeVTKCode; }

	size_t gaussePoints() const { return NodeGPCount; }

	const std::vector<DenseMatrix>& dN() const { return Node::_dN; }
	const std::vector<DenseMatrix>& N() const { return Node::_N; }
	const std::vector<double>& weighFactor() const { return Node::_weighFactor; }

	const std::vector<Property>& elementDOFs() const { return Node::_DOFElement; }
	const std::vector<Property>& faceDOFs() const { return Node::_DOFFace; }
	const std::vector<Property>& edgeDOFs() const { return Node::_DOFEdge; }
	const std::vector<Property>& pointDOFs() const { return Node::_DOFPoint; }
	const std::vector<Property>& midPointDOFs() const { return Node::_DOFMidPoint; }

private:
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
