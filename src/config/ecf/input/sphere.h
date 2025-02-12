
#ifndef SRC_CONFIG_ECF_INPUT_SPHERE_H_
#define SRC_CONFIG_ECF_INPUT_SPHERE_H_

#include "block.h"

#include <map>

namespace espreso {

struct SphereGeneratorConfiguration: public BlockGeneratorConfiguration {

    double inner_radius, outer_radius;
    size_t clusters, layers;

    std::map<std::string, std::string> nodes, edges, faces, elements;

    SphereGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_SPHERE_H_ */
