
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFROOT_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFROOT_H_

#include "xdmfelement.h"

#include <string>

namespace espreso {

class LightData: public XDMFElement {
public:
    int version;
    std::string filenane, dir;

    LightData(const std::string &filename);
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFROOT_H_ */
