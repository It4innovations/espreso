
#ifndef SRC_INPUT_OPENFOAM_PARSER_SETS_H_
#define SRC_INPUT_OPENFOAM_PARSER_SETS_H_

#include <vector>

#include "parser.h"

namespace espreso {

#define MAX_NAME_SIZE 80

struct OpenFOAMSet {
	enum class SetType {
		CELL_SET,
		FACE_SET,
		POINT_SET
	};

	char name[MAX_NAME_SIZE];
	SetType type;

	OpenFOAMSet();
	OpenFOAMSet(const std::string &name, SetType type);
};

struct OpenFOAMSets: public OpenFOAMCollectiveParser {

	OpenFOAMSets(const char *begin, const char *end): OpenFOAMCollectiveParser(begin, end) {}

	static void inspect(const std::string &path, std::vector<OpenFOAMSet> &sets);

	bool readData(OpenFOAMSet &set, std::vector<esint> &indices);
};

}



#endif /* SRC_INPUT_OPENFOAM_PARSER_SETS_H_ */
