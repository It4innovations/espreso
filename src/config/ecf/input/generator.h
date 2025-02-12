
#ifndef SRC_CONFIG_ECF_INPUT_GENERATOR_H_
#define SRC_CONFIG_ECF_INPUT_GENERATOR_H_

#include "config/description.h"

#include "grid.h"
#include "gridset.h"
#include "gridtower.h"
#include "sphere.h"

namespace espreso {

enum class INPUT_GENERATOR_SHAPE {
    GRID,
    GRID_SET,
    GRID_TOWER,
    SPHERE
};

struct InputGeneratorConfiguration: public ECFDescription {

    INPUT_GENERATOR_SHAPE shape;

    bool uniform_clusters, uniform_domains;

    GridGeneratorConfiguration grid;
    GridSetGeneratorConfiguration grid_set;
    GridTowerGeneratorConfiguration grid_tower;
    SphereGeneratorConfiguration sphere;

    InputGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_GENERATOR_H_ */
