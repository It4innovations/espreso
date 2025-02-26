
#ifndef SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_
#define SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_

#include "coordinatesystem.h"
#include "elasticityproperties.h"
#include "thermalexpansion.h"
#include "thermalconductivity.h"

namespace espreso {

struct ECFParameter;

struct MaterialBaseConfiguration: public ECFDescription {

	enum PHYSICAL_MODEL {
		THERMAL              = 1 << 0,
		STRUCTURAL_MECHANICS = 1 << 1,
		ACOUSTICS            = 1 << 2
	};

	PHYSICAL_MODEL physical_model;

	CoordinateSystemConfiguration coordinate_system;

	ECFExpression density;
	ECFExpression speed_of_sound;
	ECFExpression heat_capacity;
	ElasticityPropertiesConfiguration elasticity_properties;
	ThermalExpansionConfiguration thermal_expansion;
	ThermalConductivityConfiguration thermal_conductivity;

	MaterialBaseConfiguration(PHYSICAL_MODEL physicalModel, bool *phase_change);

protected:
	bool *_phase_change;
};

struct MaterialConfiguration: public MaterialBaseConfiguration {

	std::string name;
	std::string description;

	bool phase_change;
	size_t smooth_step_order;
	double latent_heat, transition_interval, phase_change_temperature;

	std::map<size_t, MaterialBaseConfiguration> phases;

	MaterialConfiguration(PHYSICAL_MODEL physicalModel = static_cast<PHYSICAL_MODEL>(~0));
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_ */
