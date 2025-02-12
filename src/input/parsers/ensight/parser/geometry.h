
#ifndef SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_
#define SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_

#include "basis/io/inputfile.h"

#include <string>
#include <vector>

namespace espreso {

struct MeshBuilder;

class EnsightGeometry {
    enum class Format {
        BINARY,
        ASCII,
        UNKNOWN
    };
    enum class IDs {
        OFF,
        GIVEN,
        ASSIGN,
        INGNORE,
        UNKNOWN
    };

    struct Header {
        Format format;
        IDs nodeIDs, elementIDs;
    };

    struct Coordinates {
        Coordinates(): offset(0), nn(0) {}
        Coordinates(size_t offset, int nn): offset(offset), nn(nn) {}

        size_t offset;
        int nn;
    };

    struct Elements {
        enum class Type {
            POINT,
            BAR2, BAR3,
            TRIA3, TRIA6, QUAD4, QUAD8,
            TETRA4, TETRA10, PYRAMID5, PYRAMID13, PENTA6, PENTA15, HEXA8, HEXA20,
            NSIDED, NFACED
        };

        Elements(): type(Type::POINT), offset(0), ne(0) {}
        Elements(Type type, size_t offset, int ne): type(type), offset(offset), ne(ne) {}

        Type type;
        size_t offset;
        int ne;
    };

public:
    EnsightGeometry(InputFilePack &geofile);

    void scan();
    void parse(MeshBuilder &mesh);

protected:
    void header();
    void scanBinary();
    void scanASCII();
    void parseBinary(MeshBuilder &mesh);
    void parseASCII(MeshBuilder &mesh);

    InputFilePack &_geofile;

    Header _header;
    std::vector<char> _parts;
    std::vector<Coordinates> _coordinates;
    std::vector<Elements> _elements;
};
}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_ */
