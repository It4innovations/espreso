
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_TIME_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_TIME_H_

#include "parallelFoamFile.h"
#include "input/input.h"

#include <cstddef>
#include <fstream>

namespace espreso {

struct OpenFOAMTime {

	static double value(const std::string &name);
};

}




#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_TIME_H_ */
