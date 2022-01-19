
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_

#include "foamFile.h"
#include "input/input.h"

#include <cstddef>

namespace espreso {

struct OpenFOAMVectorField: FoamFile, RawFoamFile {

	OpenFOAMVectorField(InputFile *input): RawFoamFile(input) {}

	void scan()
	{
		scanFile(*this);
	}

	void parse(ivector<_Point<esfloat> > &coordinates);
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_ */
