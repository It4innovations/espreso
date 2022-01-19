
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_

#include "foamFile.h"

#include <cstddef>

namespace espreso {

class OpenFOAMList: FoamFile, RawFoamFile {

	OpenFOAMList(InputFile *input): RawFoamFile(input) {}

	void scan()
	{
		scanFile(*this);
	}

	void parse(ivector<esint> &list);
};

}



#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_ */
