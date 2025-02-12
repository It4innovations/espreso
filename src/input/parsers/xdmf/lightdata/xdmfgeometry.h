
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFGEOMETRY_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFGEOMETRY_H_

#include "xdmfelement.h"
#include <string>

namespace espreso {

class XDMFGeometry: public XDMFElement {
public:
    enum class Type { XYZ, XY, X_Y_Z, VxVyVzm, Origin_DxDyDz, Origin_DxDy };

    std::string name;
    std::string reference;
    Type type;

    XDMFGeometry();
    void parse(XML::Element *e);
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFGEOMETRY_H_ */
