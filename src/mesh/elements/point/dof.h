
#ifndef SRC_MESH_ELEMENTS_POINT_DOF_H_
#define SRC_MESH_ELEMENTS_POINT_DOF_H_

#include "../element.h"

#define DOFNodesCount 1
#define DOFEdgeCount 0
#define DOFFacesCount 0
#define DOFGPCount 0
#define DOFCommonNodes 1
#define DOFVTKCode -1

namespace espreso {

class DOF: public Element
{

public:

	DOF(eslocal index): _index(index) {};
	Element* copy() const { return new DOF(*this); }

	eslocal nCommon() const { return DOFCommonNodes; }
	eslocal vtkCode() const { return DOFVTKCode; }
	eslocal param(Params param) const { ESINFO(GLOBAL_ERROR) << "Node has no params"; return 0; }
	void setParam(Params param, eslocal value) { ESINFO(GLOBAL_ERROR) << "Node has no params"; }
	size_t params() const { return 0; }

	size_t faces() const { return DOFFacesCount; }
	size_t edges() const { return DOFEdgeCount; }
	size_t nodes() const { return DOFNodesCount; }
	size_t coarseNodes() const { return DOFNodesCount; }
	size_t gaussePoints() const { return DOFGPCount; }

	virtual Point faceNormal(const Element *face) { ESINFO(GLOBAL_ERROR) << "Node has no face"; return Point(); }
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) { ESINFO(GLOBAL_ERROR) << "Node has no edge"; return Point(); }
	virtual Element* face(size_t index) const { ESINFO(GLOBAL_ERROR) << "Node has no face"; return NULL; }
	virtual Element* edge(size_t index) const { ESINFO(GLOBAL_ERROR) << "Node has no edge"; return NULL; }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const { ESINFO(GLOBAL_ERROR) << "Node has no neighbours"; return {}; }
	eslocal* indices() { return &_index; }
	const eslocal* indices() const { return &_index; }

	void setFace(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Node has no face"; }
	void setEdge(size_t index, Element* edge) { ESINFO(GLOBAL_ERROR) << "Node has no edge"; }
	void setFace(Element* face) { ESINFO(GLOBAL_ERROR) << "Node has no face"; }
	void setEdge(Element* edge) { ESINFO(GLOBAL_ERROR) << "Node has no edge"; }

	void fillFaces() {};
	void fillEdges() {};

private:
	eslocal _index;
};

}



#endif /* SRC_MESH_ELEMENTS_POINT_DOF_H_ */
