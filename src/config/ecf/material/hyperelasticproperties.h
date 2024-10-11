
#ifndef SRC_CONFIG_ECF_MATERIAL_HYPERELASTICPROPERTIES_H_
#define SRC_CONFIG_ECF_MATERIAL_HYPERELASTICPROPERTIES_H_

#include "tensor.h"
#include "coordinatesystem.h" // DIMENSION

namespace espreso {

struct HyperElasticPropertiesConfiguration: public ECFDescription {

	enum class MODEL {
	    KIRCHHOFF,
		NEO_HOOKEN_CMP,
		NEO_HOOKEN_INC,
		MOONEY_RIVLIN_2PARAMS,
		MOONEY_RIVLIN_3PARAMS,
		MOONEY_RIVLIN_5PARAMS,
		MOONEY_RIVLIN_9PARAMS,
		ARRUDA_BOYCE,
		BLATZ_KO_FOAM,
		GENT,
		OGDEN_1,
		OGDEN_2,
		OGDEN_3
	};

	MODEL model;

	ECFExpression E, mu;
	ECFExpression d, G, lambda;
	ECFExpression C10, C01, C11, C02, C20, C30, C21, C12, C03;

	HyperElasticPropertiesConfiguration();
};

}




#endif /* SRC_CONFIG_ECF_MATERIAL_HYPERELASTICPROPERTIES_H_ */
