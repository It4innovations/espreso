
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFDOMAIN_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFDOMAIN_H_

#include "xdmfelement.h"
#include <string>

namespace espreso {

class XDMFDomain: public XDMFElement {
public:
	std::string name;
	std::string reference;

	XDMFDomain();
	void parse(XML::Element *e);
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFDOMAIN_H_ */
