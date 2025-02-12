
#ifndef SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDSETTINGS_H_
#define SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDSETTINGS_H_

#include "block.h"
#include "blocksettings.h"

namespace espreso {

struct GridGeneratorConfiguration;
struct SphereGeneratorConfiguration;

struct GridSettings: public BlockSettings {

    GridSettings(const GridGeneratorConfiguration &configuration);
    GridSettings(const SphereGeneratorConfiguration &configuration);

    Triple<size_t> blocks, clusters;

    std::vector<bool> nonempty;
    esint body;
};

}


#endif /* SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDSETTINGS_H_ */
