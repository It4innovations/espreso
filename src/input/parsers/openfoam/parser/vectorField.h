
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_

#include "foamFile.h"
#include "input/input.h"

#include <cstddef>
#include <fstream>

namespace espreso {

struct OpenFOAMVectorField: RawFoamFile, FoamFile {

	OpenFOAMVectorField(InputFile *input): RawFoamFile(input) {}

	void scan()
	{
		scanFile(*this);
	}

	void parse(ivector<_Point<esfloat> > &coordinates);

	static esint load(const std::string &file, ivector<_Point<esfloat> > &coordinates);
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_ */
