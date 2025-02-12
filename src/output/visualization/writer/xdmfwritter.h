
#ifndef SRC_OUTPUT_VISUALIZATION_WRITER_XDMFWRITTER_H_
#define SRC_OUTPUT_VISUALIZATION_WRITER_XDMFWRITTER_H_

#include "mesh/element.h"

namespace espreso {

struct XDMFWritter {

    static int ecode(const Element::CODE &code)
    {
        switch (code) {
        case Element::CODE::POINT1:
            return 1;
        case Element::CODE::LINE2:
            return 2;
        case Element::CODE::LINE3:
            return 34;
        case Element::CODE::SQUARE4:
            return 5;
        case Element::CODE::SQUARE8:
            return 37;
        case Element::CODE::TRIANGLE3:
            return 4;
        case Element::CODE::TRIANGLE6:
            return 36;
        case Element::CODE::TETRA4:
            return 6;
        case Element::CODE::TETRA10:
            return 38;
        case Element::CODE::PYRAMID5:
            return 7;
        case Element::CODE::PYRAMID13:
            return 39;
        case Element::CODE::PRISMA6:
            return 8;
        case Element::CODE::PRISMA15:
            return 40;
        case Element::CODE::HEXA8:
            return 9;
        case Element::CODE::HEXA20:
            return 48;
        default:
            return -1;
        }
    }

    static std::string etype(const Element::CODE &code)
    {
        switch (code) {
        case Element::CODE::POINT1:
            return "Polyvertex";
        case Element::CODE::LINE2:
            return "Polyline";
        case Element::CODE::LINE3:
            return "Edge_3";
        case Element::CODE::SQUARE4:
            return "Quadrilateral";
        case Element::CODE::SQUARE8:
            return "Quad_8";
        case Element::CODE::TRIANGLE3:
            return "Triangle";
        case Element::CODE::TRIANGLE6:
            return "Tri_6";
        case Element::CODE::TETRA4:
            return "Tetrahedron";
        case Element::CODE::TETRA10:
            return "Tet_10";
        case Element::CODE::PYRAMID5:
            return "Pyramid";
        case Element::CODE::PYRAMID13:
            return "Pyramid_13";
        case Element::CODE::PRISMA6:
            return "Wedge";
        case Element::CODE::PRISMA15:
            return "Wedge_15";
        case Element::CODE::HEXA8:
            return "Hexahedron";
        case Element::CODE::HEXA20:
            return "Hex_20";
        default:
            return "Mixed";
        }
    }
};
}

#endif /* SRC_OUTPUT_VISUALIZATION_WRITER_XDMFWRITTER_H_ */
