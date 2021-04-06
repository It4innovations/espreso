
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFGRID_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFGRID_H_

#include "xdmfelement.h"
#include <string>

namespace espreso {

class XDMFGrid: public XDMFElement {
public:
	enum class Type           { Uniform, Collection, Tree, Subset };
	enum class CollectionType { Spatial, Temporal };
	enum class Section        { DataItem, All };

	std::string name;
	std::string reference;
	Type type;
	CollectionType collectiontype;
	Section section;

	XDMFGrid();
	void parse(XML::Element *e);
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFGRID_H_ */
