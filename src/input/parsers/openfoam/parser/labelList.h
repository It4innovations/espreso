
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LABELLIST_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LABELLIST_H_

#include "parallelFoamFile.h"
#include <cstddef>

namespace espreso {

struct OpenFOAMLabelList: ParallelFoamFile, RawParallelFoamFile {

	OpenFOAMLabelList(InputFile *input): RawParallelFoamFile(input) {}

	void scan()
	{
		scanFile(*this);
	}

	void parse(ivector<esint> &list);

	static FoamFileHeader load(const std::string &file, ivector<esint> &list, esint offset);
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_LABELLIST_H_ */
