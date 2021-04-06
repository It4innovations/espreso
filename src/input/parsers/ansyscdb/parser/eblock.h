
#ifndef SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_
#define SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_

#include "parser.h"

namespace espreso {

struct AnsysCDBData;
struct ET;

struct EBlock: public WorkbenchParser {
	esint NUM_NODES, Solkey, NDMAX, NDSEL;

	EBlock();
	EBlock& parse(const char* begin);

	bool readData(const std::vector<ET> &et, AnsysCDBData &mesh);

protected:
	esint lineEndSize;
	esint valueSize, valueLength;

	bool solid(const std::vector<ET> &et, AnsysCDBData &mesh);
	bool boundary(const std::vector<ET> &et, AnsysCDBData &mesh);
};

}


#endif /* SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_ */
