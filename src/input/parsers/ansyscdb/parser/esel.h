
#ifndef SRC_INPUT_WORKBENCH_PARSER_ESEL_H_
#define SRC_INPUT_WORKBENCH_PARSER_ESEL_H_

#include "parser.h"

namespace espreso {

struct ESel: public WorkbenchParser {
    enum class Type: int {
        UNKNOWN,
        S,
        R,
        A,
        U,
        ALL,
        NONE,
        INVE,
        STAT
    };

    enum class Item: int {
        UNKNOWN,
        ELEM,
        ADJ,
        CENT,
        TYPE,
        ENAME,
        MAT,
        REAL,
        ESYS,
        PART,
        LIVE,
        LAYER,
        SEC,
        STRA,
        SFE,
        BFE,
        PATH,
        ETAB
    };

    enum class Comp: int {
        UNKNOWN
    };

    Type type;
    Item item;
    Comp comp;
    esint VMIN, VMAX, VINC;
    bool KABS;

    ESel();
    ESel& parse(const char* begin);
};

}


#endif /* SRC_INPUT_WORKBENCH_PARSER_ESEL_H_ */
