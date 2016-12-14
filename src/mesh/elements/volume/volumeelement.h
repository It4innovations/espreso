
#ifndef SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_
#define SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_

#include "../element.h"

namespace espreso {

class VolumeElement: public Element
{

	friend class Mesh;

public:
	Type type() const { return Type::VOLUME; }

	virtual eslocal param(Params param) const { return _params[param]; };
	virtual void setParam(Params param, eslocal value) { _params[param] = value; }
	virtual size_t params() const { return PARAMS_SIZE; }

	virtual Element* face(size_t index) const { return _faces[index]; }
	virtual Element* edge(size_t index) const { return _edges[index]; }

	virtual void addFace(Element* face)
	{
		_faces.push_back(face);
		face->parentElements().push_back(this);
	}
	virtual void addEdge(Element* edge)
	{
		_edges.push_back(edge);
		edge->parentElements().push_back(this);
	}

protected:
	void setFace(size_t index, Element* face) { _faces[index] = face; }
	void setEdge(size_t index, Element* edge) { _edges[index] = edge; }

	eslocal _params[PARAMS_SIZE];
	std::vector<Element*> _edges;
	std::vector<Element*> _faces;
};

}

#endif /* SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_ */
