
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_

#include "config/ecf/material/material.h"

namespace espreso {

struct PhysicsConfiguration: public ECFDescription {

    enum class TYPE {
//        THERMO_ELASTICITY,
        HEAT_TRANSFER,
        STRUCTURAL_MECHANICS,
//        ACOUSTICS,
//        SHALLOW_WATER
    };

    enum class INTERPOLATION {
        LINEAR,
        QUADRATIC
    };

    enum class DISCRETIZATION {
        FEM_LOADED,
        FEM_LINEAR,
        FEM_QUADRATIC,
        FEM_TDNNS,
        BEM
    };

    int load_steps;

    // TODO: case insensitive compare
    INTERPOLATION interpolation;
    std::map<std::string, DISCRETIZATION> discretization;
    MaterialConfiguration::PHYSICAL_MODEL physical_model;

    std::map<std::string, MaterialConfiguration> materials;
    std::map<std::string, std::string> material_set;

    std::map<std::string, ECFExpression> initial_temperature, thickness;

    bool reassembling_optimization;

    PhysicsConfiguration(MaterialConfiguration::PHYSICAL_MODEL physicalModel);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_ */
