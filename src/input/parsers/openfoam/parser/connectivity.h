
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONNECTIVITY_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONNECTIVITY_H_

#include "input/input.h"

#include "containers/labelList.h"

namespace espreso {

struct OpenFOAMConnectivity {

	OpenFOAMConnectivity(InputFile *owner, InputFile *neighbors)
	: ownerFile(owner), neighborsFile(neighbors)
	{

	}

	void scan();
	void parse(OrderedUniqueFacesRegions *connectivity);

	OpenFOAMLabelList owner, neighbors;

	InputFile *ownerFile, *neighborsFile;
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONNECTIVITY_H_ */
