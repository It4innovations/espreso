
#ifndef SRC_CONFIG_ECF_MATERIAL_THERMALCONDUCTIVITY_H_
#define SRC_CONFIG_ECF_MATERIAL_THERMALCONDUCTIVITY_H_

#include "tensor.h"
#include "coordinatesystem.h" // DIMENSION

namespace espreso {

struct ThermalConductivityConfiguration: public ECFDescription {

	enum class MODEL {
		ISOTROPIC,
		DIAGONAL,
		SYMMETRIC,
		ANISOTROPIC,
	};

	MODEL model;
	DIMENSION *dimension;

	TensorConfiguration values;

	ThermalConductivityConfiguration(DIMENSION *D);

	bool needCoordinates() const;
	bool needTemperature() const;
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_THERMALCONDUCTIVITY_H_ */
