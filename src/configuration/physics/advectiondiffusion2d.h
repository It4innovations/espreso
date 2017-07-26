
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
		{ "ISOTROPIC"  , MATERIAL_MODEL::ISOTROPIC  , { "KXX" }, "Isotropic." },
		{ "DIAGONAL"   , MATERIAL_MODEL::DIAGONAL   , { "KXX", "KYY" }, "Diagonal." },
		{ "SYMMETRIC"  , MATERIAL_MODEL::SYMMETRIC  , { "KXX", "KYY", "KXY" }, "Symmetric." },
		{ "ANISOTROPIC", MATERIAL_MODEL::ANISOTROPIC, { "KXX", "KYY", "KXY", "KYX" }, "Anisotropic." }
	}));

	SUBCONFIG(CoordinateSystem, coordinate_system, "Element coordinate system.");
};

struct AdvectionDiffusion2DConfiguration: public AdvectionDiffusionConfiguration {

	SUBMAPTOMAP(size_t, std::string, std::string, thickness, "Thickness", "TIME_STEP", "Thickness settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOCONFIG(std::string, AdvectionDiffusion2DMaterial, materials, "Material description.", "MATERIAL", "Material description");

};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION2D_H_ */
