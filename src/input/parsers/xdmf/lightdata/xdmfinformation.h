
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFINFORMATION_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFINFORMATION_H_

#include "xdmfelement.h"
#include <string>

namespace espreso {

class XDMFInformation: public XDMFElement {
public:
    std::string name;
    std::string reference;
    std::string value;

    XDMFInformation();
    void parse(XML::Element *e);
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFINFORMATION_H_ */
