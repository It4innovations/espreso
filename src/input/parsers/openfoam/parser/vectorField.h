
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_

#include "parallelFoamFile.h"
#include "input/input.h"

#include <cstddef>
#include <fstream>

namespace espreso {

struct OpenFOAMVectorField: RawParallelFoamFile, ParallelFoamFile {

	OpenFOAMVectorField(InputFile *input): RawParallelFoamFile(input) {}

	void scan()
	{
		scanFile(*this);
	}

	void parse(ivector<_Point<esfloat> > &coordinates);

	static esint load(const std::string &file, ivector<_Point<esfloat> > &coordinates);
	static void load(InputFile *input, std::vector<esfloat> &data);
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_VECTORFIELD_H_ */
