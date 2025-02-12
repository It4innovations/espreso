
#ifndef SRC_INPUT_ABAQUS_PARSER_MATERIAL_H_
#define SRC_INPUT_ABAQUS_PARSER_MATERIAL_H_

#include "parser.h"
#include <map>

namespace espreso {

struct AbaqusMaterial: public AbaqusParser {
    static size_t size;
    static const char* upper;
    static const char* lower;
    static const char* sentence;

    char Name[MAX_NAME_SIZE];
    double youngs_modulus;
    double poisson_ratio;


    AbaqusMaterial();
    AbaqusMaterial& parse(const char* begin);


};

}


#endif /* SRC_INPUT_ABAQUS_COMMAND_H_ */
