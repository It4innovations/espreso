
#ifndef SRC_INPUT_OPENFOAM_PARSER_ZONES_H_
#define SRC_INPUT_OPENFOAM_PARSER_ZONES_H_

#include <vector>

#include "parser.h"

namespace espreso {

struct OpenFOAMData;
struct InputFile;

struct OpenFOAMZones: public OpenFOAMCollectiveParser {

    OpenFOAMZones(InputFile &pfile);

    int getZones();
    void synchronize(int zones, std::vector<char> &names, std::vector<size_t> &offsets);
    void readData(std::vector<esint> &indices, size_t begin, size_t end);

    bool readPoints(OpenFOAMData &data);
    bool readFaces(OpenFOAMData &data);
    bool readCells(OpenFOAMData &data);

    InputFile &_pfile;
};

}


#endif /* SRC_INPUT_OPENFOAM_PARSER_ZONES_H_ */
