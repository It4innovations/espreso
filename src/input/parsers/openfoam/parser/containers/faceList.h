
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FACELIST_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FACELIST_H_

#include "foamFile.h"
#include "mesh/element.h"

#include <cstddef>

namespace espreso {

struct OpenFOAMFaceList: FoamFile, FoamFileDistribution {

	void scan(InputFile *file)
	{
		FoamFile::scan(file, *this);
	}

	void parse(ivector<Element::CODE> &type, ivector<esint> &enodes);
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FACELIST_H_ */
