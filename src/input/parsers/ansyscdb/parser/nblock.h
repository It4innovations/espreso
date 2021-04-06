
#ifndef SRC_INPUT_WORKBENCH_COMMAND_H_
#define SRC_INPUT_WORKBENCH_COMMAND_H_

#include "parser.h"
#include "basis/containers/point.h"

namespace espreso {

struct AnsysCDBData;

struct NBlock: public WorkbenchParser {
	esint NUMFIELD, Solkey, NDMAX, NDSEL;

	esint lineSize, lineEndSize;
	esint indexSize, indexLength, valueSize, valueLength;

	NBlock();
	NBlock& parse(const char* begin);
	bool readData(AnsysCDBData &mesh);
};

}


#endif /* SRC_INPUT_WORKBENCH_COMMAND_H_ */
