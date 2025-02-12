
#ifndef SRC_INPUT_WORKBENCH_WORKBENCH_H_
#define SRC_INPUT_WORKBENCH_WORKBENCH_H_

#include "input/meshbuilder.h"
#include "basis/io/inputfile.h"

#include <cstddef>
#include <vector>

namespace espreso {

struct InputConfiguration;

struct NBlock;
struct EBlock;
struct CMBlock;
struct ET;
struct ESel;
struct NSel;
struct CM;
struct BlockEnd;

struct AnsysCDBData: public MeshBuilder {
    std::vector<int> et;
};

class AnsysCDBLoader: public AnsysCDBData {
public:
    AnsysCDBLoader(const InputConfiguration &configuration);
    void load();

protected:
    void scan();
    void parse();

    const InputConfiguration &_configuration;

    std::vector<NBlock> _NBlocks;
    std::vector<EBlock> _EBlocks;
    std::vector<CMBlock> _CMBlocks;
    std::vector<ET> _ET;
    std::vector<ESel> _ESel;
    std::vector<NSel> _NSel;
    std::vector<CM> _CM;
    std::vector<BlockEnd> _blockEnds;

    InputFilePack _file;
};

}



#endif /* SRC_INPUT_WORKBENCH_WORKBENCH_H_ */
