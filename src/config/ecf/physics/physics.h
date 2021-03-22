
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_

#include "config/ecf/material/material.h"

namespace espreso {

struct PhysicsConfiguration: public ECFDescription {

	enum class TYPE {
		THERMO_ELASTICITY_2D,
		THERMO_ELASTICITY_3D,
		HEAT_TRANSFER_2D,
		HEAT_TRANSFER_3D,
		STRUCTURAL_MECHANICS_2D,
		STRUCTURAL_MECHANICS_3D,
		ACOUSTICS_2D,
		ACOUSTICS_3D,
		SHALLOW_WATER_2D
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
	DIMENSION dimension;
	MaterialConfiguration::PHYSICAL_MODEL physical_model;

	std::map<std::string, MaterialConfiguration> materials;
	std::map<std::string, std::string> material_set;

	std::map<std::string, ECFExpression> initial_temperature, thickness;

	bool contact_interfaces;

	PhysicsConfiguration(DIMENSION dim, MaterialConfiguration::PHYSICAL_MODEL physicalModel);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_ */
