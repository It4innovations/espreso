
#ifndef SRC_MESH_ELEMENTS_POINT_NODE_H_
#define SRC_MESH_ELEMENTS_POINT_NODE_H_

#include "../element.h"

#define NodeNodesCount 0
#define NodeFacesCount 0
#define NodeGPCount 0
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

	Node(): Element(_nparams) { };

	Element* copy() const
	{
		return new Node(*this);
	}

	eslocal vtkCode() const
	{
		return NodeVTKCode;
	}

	const eslocal* indices() const
	{
		return indices();
	}

	size_t size() const
	{
		return NodeNodesCount;
	}

	size_t coarseSize() const
	{
		return NodeNodesCount;
	}

	size_t gpSize() const
	{
		return NodeGPCount;
	}

	size_t faces() const
	{
		return NodeFacesCount;
	}

	const std::vector<DenseMatrix>& dN() const
	{
		return Node::_dN;
	}

	const std::vector<DenseMatrix>&  N() const
	{
		return Node::_N;
	}

	const std::vector<double>& weighFactor() const
	{
		return Node::_weighFactor;
	}

	const std::vector<Property>& DOFElement() const
	{
		return Node::_DOFElement;
	}

	const std::vector<Property>& DOFFace() const
	{
		return Node::_DOFFace;
	}

	const std::vector<Property>& DOFEdge() const
	{
		return Node::_DOFEdge;
	}

	const std::vector<Property>& DOFPoint() const
	{
		return Node::_DOFPoint;
	}

	const std::vector<Property>& DOFMidPoint() const
	{
		return Node::_DOFMidPoint;
	}

	eslocal nCommon() const
	{
		return 1;
	}

	std::vector<eslocal> getNeighbours(size_t nodeIndex) const
	{
		ESINFO(GLOBAL_ERROR) << "Node does not have neighbours.";
		return std::vector<eslocal>();
	}

	std::vector<eslocal> getFace(size_t face) const
	{
		ESINFO(GLOBAL_ERROR) << "Node does not have faces.";
		return std::vector<eslocal>();
	}

	Element* getFullFace(size_t face) const
	{
		ESINFO(GLOBAL_ERROR) << "Node does not have faces.";
		return NULL;
	}

	Element* getCoarseFace(size_t face) const
	{
		ESINFO(GLOBAL_ERROR) << "Node does not have faces.";
		return NULL;
	}

protected:

	eslocal* indices()
	{
		ESINFO(GLOBAL_ERROR) << "Node does not have indices.";
		return NULL;
	}

private:
	static size_t _counter;

	static eslocal _nparams[PARAMS_SIZE];

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
