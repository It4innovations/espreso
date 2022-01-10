
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_

#include "foamFile.h"

#include <cstddef>

namespace espreso {

class OpenFOAMList: FoamFile, FoamFileDistribution {

	void scan(InputFile *file)
	{
		FoamFile::scan(file, *this);
	}

	void parse(ivector<esint> &list);
};

}



#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LIST_H_ */
