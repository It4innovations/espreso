
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FACELIST_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FACELIST_H_

#include "parallelFoamFile.h"
#include "mesh/element.h"

#include <cstddef>

namespace espreso {

struct OpenFOAMFaceList: ParallelFoamFile, RawParallelFoamFile {

	OpenFOAMFaceList(InputFile *input): RawParallelFoamFile(input) {}

	void scan()
	{
		scanFile(*this);
	}

	void parse(ivector<Element::CODE> &type, ivector<esint> &enodes);

	static FoamFileHeader load(const std::string &file, ivector<Element::CODE> &type, ivector<esint> &enodes, esint offset);
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FACELIST_H_ */
