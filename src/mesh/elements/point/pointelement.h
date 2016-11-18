
#ifndef SRC_MESH_ELEMENTS_POINT_POINTELEMENT_H_
#define SRC_MESH_ELEMENTS_POINT_POINTELEMENT_H_

#include "../element.h"

#define PointNodesCount 1
#define PointEdgeCount 0
#define PointFacesCount 0
#define PointCommonNodes 1

namespace espreso {

class PointElement: public Element
{

public:
	PointElement(eslocal index): _index(index) {};

	Type type() const { return Type::POINT; }

	eslocal param(Params param) const { ESINFO(GLOBAL_ERROR) << "Point element has no params"; return 0; }
	void setParam(Params param, eslocal value) { ESINFO(GLOBAL_ERROR) << "Point element has no params"; }
	size_t params() const { return 0; }

	eslocal nCommon() const { return PointCommonNodes; }
	size_t faces() const { return PointFacesCount; }
	size_t edges() const { return PointEdgeCount; }
	size_t nodes() const { return PointNodesCount; }
	size_t coarseNodes() const { return PointNodesCount; }

	virtual Element* face(size_t index) const
	{
		ESINFO(GLOBAL_ERROR) << "Point element has no face";
		return NULL;
	}
	virtual Element* edge(size_t index) const
	{
		ESINFO(GLOBAL_ERROR) << "Point element has no edge";
		return NULL;
	}

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const
	{
		ESINFO(GLOBAL_ERROR) << "Point element has no neighbours";
		return {};
	}

	void setFace(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Point element has no face"; }
	void setEdge(size_t index, Element* edge) { ESINFO(GLOBAL_ERROR) << "Point element has no edge"; }
	void setFace(Element* face) { ESINFO(GLOBAL_ERROR) << "Point element has no face"; }
	void setEdge(Element* edge) { ESINFO(GLOBAL_ERROR) << "Point element has no edge"; }

	void fillFaces() { ESINFO(GLOBAL_ERROR) << "Call fill faces on Point element"; }
	void fillEdges() { ESINFO(GLOBAL_ERROR) << "Call fill edges on Point element"; }

	eslocal* indices() { return &_index; }
	const eslocal* indices() const { return &_index; }

private:
	eslocal _index;
};

}

#endif /* SRC_MESH_ELEMENTS_POINT_POINTELEMENT_H_ */
