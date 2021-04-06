
#ifndef SRC_INPUT_WORKBENCH_PARSER_CMBLOCK_H_
#define SRC_INPUT_WORKBENCH_PARSER_CMBLOCK_H_

#include "parser.h"

namespace espreso {

struct CMBlock: public WorkbenchParser {
	enum class Entity: int {
		UNKNOWN,
		NODE,
		ELEMENT,
		KP,
		LINE,
		AREA,
		VOLU
	};

	char name[MAX_NAME_SIZE];
	Entity entity;
	esint NUMITEMS;

	esint lineSize, lineEndSize;
	esint valueSize, valueLength;

	CMBlock();
	CMBlock& parse(const char* begin);

	bool readData(std::vector<esint> &indices);
};

}


#endif /* SRC_INPUT_WORKBENCH_PARSER_CMBLOCK_H_ */
