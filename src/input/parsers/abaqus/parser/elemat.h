#ifndef SRC_INPUT_ABAQUS_PARSER_ELEMAT_H_
#define SRC_INPUT_ABAQUS_PARSER_ELEMAT_H_

#include "parser.h"

#include <map>
#include <string>

namespace espreso {

struct SSection;
struct AbaqusMaterial;

struct Elemat: public AbaqusParser {

    Elemat();
    Elemat& create_dict(const std::vector<SSection> &ssection, const std::vector<AbaqusMaterial> &material);

    std::map<std::string, std::string> strig11;
    std::string name_mat;
    std::string namee;
    std::string namee1;
    std::vector<std::string> elset_mat_dict;

};

}


#endif /* SRC_INPUT_ABAQUS_COMMAND_H_ */
