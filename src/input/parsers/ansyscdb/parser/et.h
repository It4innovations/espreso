
#ifndef SRC_INPUT_WORKBENCH_PARSER_ET_H_
#define SRC_INPUT_WORKBENCH_PARSER_ET_H_

#include "parser.h"

namespace espreso {

struct ET: public WorkbenchParser {
    enum class ETYPE {
        D2SOLID_4NODES,
        D2SOLID_6NODES,
        D2SOLID_8NODES,

        D3SOLID_4NODES,
        D3SOLID_8NODES,
        D3SOLID_10NODES,
        D3SOLID_20NODES,

        D2SURFACE,
        D3SURFACE,

        TARGET,
        CONTACT,
        PRETS,

        UNIVERSAL,
        SKIP,

        UNKNOWN
    };

    esint id, type;

    ET();
    ET& parse(const char* begin);

    ETYPE etype() const
    {
        switch (type) {
        case  35: return ETYPE::D2SOLID_6NODES;
        case  55: return ETYPE::D2SOLID_4NODES;
        case  77: return ETYPE::D2SOLID_8NODES;
        case 162: return ETYPE::D2SOLID_4NODES;
        case 182: return ETYPE::D2SOLID_4NODES;
        case 183: return ETYPE::D2SOLID_8NODES;

        case  45: return ETYPE::D3SOLID_8NODES;
        case  70: return ETYPE::D3SOLID_8NODES;
        case  87: return ETYPE::D3SOLID_10NODES;
        case  90: return ETYPE::D3SOLID_20NODES;
        case 164: return ETYPE::D3SOLID_8NODES;
        case 168: return ETYPE::D3SOLID_10NODES;
        case 185: return ETYPE::D3SOLID_8NODES;
        case 186: return ETYPE::D3SOLID_20NODES;
        case 187: return ETYPE::D3SOLID_10NODES;
        case 278: return ETYPE::D3SOLID_8NODES;
        case 279: return ETYPE::D3SOLID_20NODES;
        case 285: return ETYPE::D3SOLID_4NODES;

        case 151: return ETYPE::D2SURFACE;
        case 152: return ETYPE::D3SURFACE;
        case 153: return ETYPE::D2SURFACE;
        case 154: return ETYPE::D3SURFACE;
        case 155: return ETYPE::D3SURFACE;
        case 156: return ETYPE::D3SURFACE;

        case 170: return ETYPE::TARGET; // TARGET
        case 174: return ETYPE::CONTACT; // CONTACT
        case 179: return ETYPE::PRETS; // PRETS

        case 200: return ETYPE::UNIVERSAL;

        case 201: return ETYPE::SKIP;

        default: return ETYPE::UNKNOWN;
        }
    }
};

}



#endif /* SRC_INPUT_WORKBENCH_PARSER_ET_H_ */
