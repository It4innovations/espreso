
#ifndef SRC_CONFIG_ECF_MATERIAL_PLASTICITYPROPERTIES_H_
#define SRC_CONFIG_ECF_MATERIAL_PLASTICITYPROPERTIES_H_

#include "tensor.h"
#include "coordinatesystem.h" // DIMENSION

namespace espreso {

struct PlasticityPropertiesConfiguration: public ECFDescription {

	DIMENSION *dimension;

	PlasticityPropertiesConfiguration(DIMENSION *D);
};

}

#endif /* SRC_CONFIG_ECF_MATERIAL_PLASTICITYPROPERTIES_H_ */
