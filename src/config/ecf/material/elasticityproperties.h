
#ifndef SRC_CONFIG_ECF_MATERIAL_ELASTICITYPROPERTIES_H_
#define SRC_CONFIG_ECF_MATERIAL_ELASTICITYPROPERTIES_H_

#include "tensor.h"

namespace espreso {

struct ElasticityPropertiesConfiguration: public ECFDescription {

    enum class MODEL {
        ISOTROPIC,
//        ORTHOTROPIC,
//        ANISOTROPIC
    };

    enum class MATERIAL_MODEL {
        KIRCHHOFF,
        NEO_HOOKEAN_CMP,
//        NEO_HOOKEAN_INC,
//        MOONEY_RIVLIN_2PARAMS,
//        MOONEY_RIVLIN_3PARAMS,
//        MOONEY_RIVLIN_5PARAMS,
//        MOONEY_RIVLIN_9PARAMS,
//        ARRUDA_BOYCE,
//        BLATZ_KO_FOAM,
//        GENT,
//        OGDEN_1,
//        OGDEN_2,
//        OGDEN_3,
        BONET_WOOD
    };

    MODEL model;
    MATERIAL_MODEL material_model;

    TensorConfiguration poisson_ratio;
    TensorConfiguration young_modulus;
    TensorConfiguration shear_modulus;
    TensorConfiguration anisotropic;

    ECFExpression sigma;
    ECFExpression Hisotropic;

    ElasticityPropertiesConfiguration();

    bool needCoordinates() const;
    bool needTemperature() const;
};

}

#endif /* SRC_CONFIG_ECF_MATERIAL_ELASTICITYPROPERTIES_H_ */
