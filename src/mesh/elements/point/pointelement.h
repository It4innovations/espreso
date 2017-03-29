
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

	virtual size_t filledFaces() const { return 0; }
	virtual size_t filledEdges() const { return 0; }

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

	void addFace(Element* face) { ESINFO(GLOBAL_ERROR) << "Point element has no face"; }
	void addEdge(Element* edge) { ESINFO(GLOBAL_ERROR) << "Point element has no edge"; }

	Element* addFace(const std::vector<eslocal> &nodes) { return NULL; }

	const std::vector<eslocal>& faceNodes(size_t index) const
	{
		static std::vector<eslocal> _facesNodes;
		return _facesNodes;
	}

	const std::vector<eslocal>& edgeNodes(size_t index) const
	{
		static std::vector<eslocal> _edgesNodes;
		return _edgesNodes;
	}

	const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Point has no base functions for face."; exit(1); }
	const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Point volume has no base functions for face."; exit(1); }
	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Point volume has no base functions for edge."; exit(1); }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Point volume has no base functions for edge."; exit(1); }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const
	{
		ESINFO(GLOBAL_ERROR) << "Point element has no neighbours";
		return {};
	}

	void setFace(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Point element has no face"; }
	void setEdge(size_t index, Element* edge) { ESINFO(GLOBAL_ERROR) << "Point element has no edge"; }

	size_t fillFaces() { ESINFO(GLOBAL_ERROR) << "Call fill faces on Point element"; return 0; }
	size_t fillEdges() { ESINFO(GLOBAL_ERROR) << "Call fill edges on Point element"; return 0; }

	eslocal* indices() { return &_index; }
	const eslocal* indices() const { return &_index; }

private:
	eslocal _index;
};

}

#endif /* SRC_MESH_ELEMENTS_POINT_POINTELEMENT_H_ */
