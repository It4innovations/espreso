
#ifndef SRC_CONFIG_ECF_INPUT_GRIDSET_H_
#define SRC_CONFIG_ECF_INPUT_GRIDSET_H_

#include "grid.h"

namespace espreso {

struct GridSetGeneratorConfiguration: public ECFDescription {

    std::map<esint, GridGeneratorConfiguration> grids;

    GridSetGeneratorConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_INPUT_GRIDSET_H_ */
