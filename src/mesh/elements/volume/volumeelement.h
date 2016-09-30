
#ifndef SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_
#define SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_

#include "../element.h"

namespace espreso {

class VolumeElement: public Element
{

public:
	Type type() const { return Type::VOLUME; }

	virtual eslocal param(Params param) const { return _params[param]; };
	virtual void setParam(Params param, eslocal value) { _params[param] = value; }
	virtual size_t params() const { return PARAMS_SIZE; }

protected:
	eslocal _params[PARAMS_SIZE];
	std::vector<Element*> _edges;
	std::vector<Element*> _faces;
};

}

#endif /* SRC_MESH_ELEMENTS_VOLUME_VOLUMEELEMENT_H_ */
