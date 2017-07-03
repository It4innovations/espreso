
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION2D_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION2D_H_

#include "advectiondiffusion.h"

namespace espreso {

struct AdvectionDiffusion2DMaterial: public Configuration {

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::DENSITY>                , density, "Density", {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::HEAT_CAPACITY>          , Cp     , "Termal capacity."       , {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX>, KXX    , "Termal conductivity XX.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY>, KYY    , "Termal conductivity YY.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY>, KXY    , "Termal conductivity XY.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX>, KYX    , "Termal conductivity YX.", {"1"});

	OPTION(MATERIAL_MODEL, model, "Material model", MATERIAL_MODEL::ISOTROPIC, OPTIONS({
		{ "ISOTROPIC"  , MATERIAL_MODEL::ISOTROPIC  , "Isotropic." },
		{ "DIAGONAL"   , MATERIAL_MODEL::DIAGONAL   , "Diagonal." },
		{ "SYMMETRIC"  , MATERIAL_MODEL::SYMMETRIC  , "Symmetric." },
		{ "ANISOTROPIC", MATERIAL_MODEL::ANISOTROPIC, "Anisotropic." }
	}));

	SUBCONFIG(CoordinateSystem, coordinate_system, "Element coordinate system.");
};

struct AdvectionDiffusion2DConfiguration: public AdvectionDiffusionConfiguration {

	SUBMAPTOMAP(size_t, std::string, std::string, thickness, "Thickness", "1", "Thickness settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOCONFIG(std::string, AdvectionDiffusion2DMaterial, materials, "Material description.", "<MATERIAL_NAME>", "Material description");

};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION2D_H_ */
