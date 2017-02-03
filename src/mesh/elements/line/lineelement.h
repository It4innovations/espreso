
#ifndef SRC_MESH_ELEMENTS_LINE_LINEELEMENT_H_
#define SRC_MESH_ELEMENTS_LINE_LINEELEMENT_H_

#include "../element.h"

namespace espreso {

class LineElement: public Element
{

public:
	Type type() const { return Type::LINE; }

	eslocal param(Params param) const { return _params[param]; }
	void setParam(Params param, eslocal value) { _params[param] = value; }
	size_t params() const { return _params.size(); }

	virtual size_t filledFaces() const { return 0; }
	virtual size_t filledEdges() const { return 0; }

	size_t faces() const { return 0; }
	size_t edges() const { return 0; }

	virtual Element* face(size_t index) const
	{
		ESINFO(GLOBAL_ERROR) << "Line element has no face";
		return NULL;
	}
	virtual Element* edge(size_t index) const
	{
		ESINFO(GLOBAL_ERROR) << "Line element has no edge";
		return NULL;
	}

	void addFace(Element* face) { ESINFO(GLOBAL_ERROR) << "Line element has no face"; }
	void addEdge(Element* edge) { ESINFO(GLOBAL_ERROR) << "Line element has no edge"; }

protected:
	void setFace(size_t index, Element* face) { ESINFO(GLOBAL_ERROR) << "Line element has no face"; }
	void setEdge(size_t index, Element* edge) { ESINFO(GLOBAL_ERROR) << "Line element has no edge"; }

	size_t fillFaces() { ESINFO(GLOBAL_ERROR) << "Call fill faces on line element."; return 0; }
	size_t fillEdges() { ESINFO(GLOBAL_ERROR) << "Call fill edges on line element."; return 0; }

	std::vector<eslocal> _params;
};

}

#endif /* SRC_MESH_ELEMENTS_LINE_LINEELEMENT_H_ */
