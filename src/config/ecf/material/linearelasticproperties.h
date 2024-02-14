
#ifndef SRC_CONFIG_ECF_MATERIAL_LINEARELASTICPROPERTIES_H_
#define SRC_CONFIG_ECF_MATERIAL_LINEARELASTICPROPERTIES_H_

#include "tensor.h"
#include "coordinatesystem.h" // DIMENSION

namespace espreso {

struct LinearElasticPropertiesConfiguration: public ECFDescription {

	enum class MODEL {
		ISOTROPIC,
		ORTHOTROPIC,
		ANISOTROPIC
	};

	bool orientation;

	MODEL model;

	TensorConfiguration poisson_ratio;
	TensorConfiguration young_modulus;
	TensorConfiguration shear_modulus;
	TensorConfiguration anisotropic;

	LinearElasticPropertiesConfiguration();

	bool needCoordinates() const;
	bool needTemperature() const;
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_LINEARELASTICPROPERTIES_H_ */
