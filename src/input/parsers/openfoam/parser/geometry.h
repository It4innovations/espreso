
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_GEOMETRY_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_GEOMETRY_H_

#include "input/input.h"

#include "containers/faceList.h"
#include "containers/vectorField.h"

namespace espreso {

struct OpenFOAMGeometry {

	OpenFOAMGeometry(InputFile *points, InputFile *faces)
	: pointsFile(points), facesFile(faces)
	{

	}

	void scan();
	void parse(OrderedUniqueNodes *nodes, OrderedUniqueFaces *faces);

	OpenFOAMVectorField points;
	OpenFOAMFaceList faces;

	InputFile *pointsFile, *facesFile;
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_GEOMETRY_H_ */
