
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFTOPOLOGY_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFTOPOLOGY_H_

#include "xdmfelement.h"
#include <string>

namespace espreso {

class XDMFTopology: public XDMFElement {
public:
    enum class Type: int {
        Polyvertex,
        Polyline,
        Polygon,
        Triangle,
        Quadrilateral,
        Tetrahedron,
        Pyramid,
        Wedge,
        Hexahedron,
        Edge_3,
        Triangle_6,
        Quadrilateral_8,
        Tetrahedron_10,
        Pyramid_13,
        Wedge_15,
        Hexahedron_20,
        Mixed,
        SMesh2D,
        RectMesh2D,
        CoRectMesh2D,
        SMesh3D,
        RectMesh3D,
        CoRectMesh3D
    };

    std::string name;
    std::string reference;
    Type type;
    int nodeperelement;
    int numberofelement;
    int dimension;
    int order;

    XDMFTopology();
    void parse(XML::Element *e);
};

}

#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFTOPOLOGY_H_ */
