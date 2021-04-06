
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFTIME_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFTIME_H_

#include "xdmfelement.h"

namespace espreso {

class XDMFTime: public XDMFElement {
public:
	enum class Type { Single, HyperSlab, List, Range };

	Type type;
	double value;

	XDMFTime();
	void parse(XML::Element *e);
};

}



#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFTIME_H_ */
