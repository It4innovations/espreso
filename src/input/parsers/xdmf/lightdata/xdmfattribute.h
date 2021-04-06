
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFATTRIBUTE_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFATTRIBUTE_H_

#include "xdmfelement.h"
#include <string>

namespace espreso {

class XDMFAttribute: public XDMFElement {
public:
	enum class Type          { Scalar, Vector, Tensor, Tensor6, Matrix, GlobalID };
	enum class Center        { Node, Cell, Grid, Face, Edge, Other };
	enum class ItemType      {};
	enum class ElementFamily {};
	enum class ElementDegree {};
	enum class ElementCell   {};

	std::string name;
	std::string reference;
	Type type;
	Center center;
	ItemType itemtype;
	ElementFamily efaminly;
	ElementDegree edegree;
	ElementCell ecell;

	XDMFAttribute();
	void parse(XML::Element *e);
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFATTRIBUTE_H_ */
