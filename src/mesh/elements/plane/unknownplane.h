
#ifndef SRC_MESH_ELEMENTS_PLANE_UNKNOWNPLANE_H_
#define SRC_MESH_ELEMENTS_PLANE_UNKNOWNPLANE_H_

#include "planeelement.h"

namespace espreso {

class UnknownPlane: public PlaneElement
{
public:
	UnknownPlane(std::vector<eslocal> &nodes): _nodes(nodes) {};
	Element* copy() const { return new UnknownPlane(*this); }

	eslocal nCommon() const { return _nodes.size() > 4 ? 3 : 2; }
	eslocal vtkCode() const { ESINFO(GLOBAL_ERROR) << "Want VTK of unknown plane element."; return -1; }

	size_t edges() const { return _edges.size(); }
	size_t nodes() const { return _nodes.size(); }
	size_t coarseNodes() const { return _nodes.size(); }
	size_t gaussePoints() const { ESINFO(GLOBAL_ERROR) << "Unknown plane has no gausse points."; return 0; }

	virtual Point faceNormal(const Element *face) const
	{
		ESINFO(GLOBAL_ERROR) << "Call normal of unknown plane element.";
		return Point();
	}
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) const
	{
		ESINFO(GLOBAL_ERROR) << "Call normal of unknown plane element.";
		return Point();
	}

	virtual Element* edge(size_t index) const { return _edges[index]; }

	const std::vector<DenseMatrix>& dN() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<DenseMatrix>& N() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<double>& weighFactor() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }

	const std::vector<Property>& elementDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& faceDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& edgeDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& pointDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& midPointDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const
	{
		ESINFO(GLOBAL_ERROR) << "Call neighbour of unknown plane element node.";
		return std::vector<eslocal>();
	}
	eslocal* indices() { return _nodes.data(); }
	const eslocal* indices() const { return _nodes.data(); }

	void setEdge(size_t index, Element* edge) { _edges[index] = edge; }
	void setEdge(Element* edge) { _edges.push_back(edge); };

	void fillEdges();

private:
	std::vector<eslocal> &_nodes;
};


}

#endif /* SRC_MESH_ELEMENTS_PLANE_UNKNOWNPLANE_H_ */
