
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_POLYBOUNDARYMESH_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_POLYBOUNDARYMESH_H_

#include "foamfileheader.h"

#include <string>

namespace espreso {

struct OpenFOAMPolyBoundaryMesh {
	enum OpenFOAMBoundaryType {
		patch, wall, processor
	};

	struct OpenFOAMBoundary {
		char name[64];
		OpenFOAMBoundaryType type;
		esint nFaces, startFace;
		int neighbor;
	};

	OpenFOAMPolyBoundaryMesh(const std::string &file);

	std::vector<OpenFOAMBoundary> boundaries;
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_POLYBOUNDARYMESH_H_ */
