
#ifndef SRC_INPUT_OPENFOAM_PARSER_POINTS_H_
#define SRC_INPUT_OPENFOAM_PARSER_POINTS_H_

#include "parser.h"
#include "basis/containers/point.h"

#include <vector>

namespace espreso {

struct OpenFOAMPoints: public OpenFOAMCollectiveParser {

    OpenFOAMPoints(const char *begin, const char *end): OpenFOAMCollectiveParser(begin, end) {}

    bool readData(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor);
};

}


#endif /* SRC_INPUT_OPENFOAM_PARSER_POINTS_H_ */
