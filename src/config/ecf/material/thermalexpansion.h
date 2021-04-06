
#ifndef SRC_CONFIG_ECF_MATERIAL_THERMALEXPANSION_H_
#define SRC_CONFIG_ECF_MATERIAL_THERMALEXPANSION_H_

#include "tensor.h"
#include "coordinatesystem.h" // DIMENSION

namespace espreso {

struct ThermalExpansionConfiguration: public ECFDescription {

	enum class MODEL {
		ISOTROPIC,
		ORTHOTROPIC
	};

	MODEL model;
	DIMENSION *dimension;

	TensorConfiguration thermal_expansion;

	ThermalExpansionConfiguration(DIMENSION *D);
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_THERMALEXPANSION_H_ */
