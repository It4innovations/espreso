
#ifndef SRC_INPUT_WORKBENCH_PARSER_BLOCKEND_H_
#define SRC_INPUT_WORKBENCH_PARSER_BLOCKEND_H_

#include "parser.h"

namespace espreso {

struct BlockEnd: public WorkbenchParser {
    BlockEnd();
    BlockEnd& parse(const char* begin);
};

}



#endif /* SRC_INPUT_WORKBENCH_PARSER_BLOCKEND_H_ */
