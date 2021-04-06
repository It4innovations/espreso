
#ifndef SRC_INPUT_OPENFOAM_PARSER_FACES_H_
#define SRC_INPUT_OPENFOAM_PARSER_FACES_H_

#include <vector>

#include "parser.h"

namespace espreso {

struct OpenFOAMData;

struct OpenFOAMFaces: public OpenFOAMCollectiveParser {

	OpenFOAMFaces(const char *begin, const char *end): OpenFOAMCollectiveParser(begin, end) {}

	bool readFaces(OpenFOAMData &data);
	bool readParents(std::vector<esint> &data);
};

}



#endif /* SRC_INPUT_OPENFOAM_PARSER_FACES_H_ */
