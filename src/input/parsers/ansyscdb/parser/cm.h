
#ifndef SRC_INPUT_WORKBENCH_PARSER_CM_H_
#define SRC_INPUT_WORKBENCH_PARSER_CM_H_

#include "parser.h"

#include <map>

namespace espreso {

struct ESel;
struct AnsysCDBData;
struct MeshERegion;

struct NSel;
struct MeshNRegion;

struct CM: public WorkbenchParser {
    enum class Entity: int {
        VOLUME,
        AREA,
        LINE,
        KP,
        ELEMENTS,
        NODES
    };

    char name[MAX_NAME_SIZE];
    Entity entity;

    CM();
    CM& parse(const char* begin);

    bool addRegion(
            const AnsysCDBData &mesh,
            const std::vector<ESel> &esel, std::map<std::string, std::vector<esint> > &eregions,
            const std::vector<NSel> &nsel, std::map<std::string, std::vector<esint> > &nregions);

protected:
    bool addElementRegion(const AnsysCDBData &mesh, const std::vector<ESel> &esel, std::map<std::string, std::vector<esint> > &eregions);
    bool addNodeRegion(const std::vector<NSel> &nsel, std::map<std::string, std::vector<esint> > &nregions);
};

}



#endif /* SRC_INPUT_WORKBENCH_PARSER_CM_H_ */
