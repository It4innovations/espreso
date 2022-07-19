
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_

#include "parallelFoamFile.h"
#include <cstddef>

namespace espreso {

class OpenFOAMList: ParallelFoamFile, RawParallelFoamFile {

	OpenFOAMList(InputFile *input): RawParallelFoamFile(input) {}

	void scan()
	{
		scanFile(*this);
	}

	void parse(ivector<esint> &list);
};

}



#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_ */
