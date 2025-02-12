
#ifndef SRC_INPUT_ABAQUS_PARSER_ESET_H_
#define SRC_INPUT_ABAQUS_PARSER_ESET_H_

#include "parser.h"

namespace espreso {

struct Eset: public AbaqusParser {
    static size_t size;
    static const char* upper;
    static const char* lower;
    static const char* sentence;
    char NAME[MAX_NAME_SIZE];

    esint NUMFIELD, Solkey, NDMAX, NDSEL;

    esint lineSize, lineEndSize;
    esint indexSize, indexLength, valueSize, valueLength;

    Eset();
    Eset& parse(const char* begin);

    bool readData(std::vector<esint> &indices);


};

}


#endif /* SRC_INPUT_ABAQUS_COMMAND_H_ */
