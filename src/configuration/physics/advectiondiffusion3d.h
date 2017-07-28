
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION3D_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION3D_H_

#include "advectiondiffusion.h"
#include "discretization.h"

namespace espreso {

struct AdvectionDiffusion3DMaterial: public Configuration {

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::DENSITY>                , density, "Density", {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::HEAT_CAPACITY>          , Cp     , "Termal capacity."       , {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX>, KXX    , "Termal conductivity XX.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY>, KYY    , "Termal conductivity YY.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ>, KZZ    , "Termal conductivity ZZ.", {"1"});

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY>, KXY    , "Termal conductivity XY.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ>, KXZ    , "Termal conductivity XZ.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ>, KYZ    , "Termal conductivity YZ.", {"1"});

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX>, KYX    , "Termal conductivity YX.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZX>, KZX    , "Termal conductivity ZX.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZY>, KZY    , "Termal conductivity ZY.", {"1"});

	OPTION(MATERIAL_MODEL, model, "Material model", MATERIAL_MODEL::ISOTROPIC, OPTIONS({
		{ "ISOTROPIC"  , MATERIAL_MODEL::ISOTROPIC  , { "KXX" }, "Isotropic." },
		{ "DIAGONAL"   , MATERIAL_MODEL::DIAGONAL   , { "KXX", "KYY", "KZZ" }, "Diagonal." },
		{ "SYMMETRIC"  , MATERIAL_MODEL::SYMMETRIC  , { "KXX", "KYY", "KZZ", "KXY", "KXZ", "KYZ" }, "Symmetric." },
		{ "ANISOTROPIC", MATERIAL_MODEL::ANISOTROPIC, { "KXX", "KYY", "KZZ", "KXY", "KXZ", "KYZ", "KYX", "KZX", "KZY" }, "Anisotropic." }
	}));

	SUBCONFIG(CoordinateSystem, coordinate_system, "Element coordinate system.");
};

struct AdvectionDiffusion3DConfiguration: public AdvectionDiffusionConfiguration {

	SUBMAPTOCONFIG(std::string, AdvectionDiffusion3DMaterial, materials, "Material description.", "MATERIAL", "Material description");

	OPTION(DISCRETIZATION, discretization, "Discretization type", DISCRETIZATION::FEM, OPTIONS({
		{ "FEM", DISCRETIZATION::FEM, "Finite element method." },
		{ "BEM", DISCRETIZATION::BEM, "Boundary element method." }
	}));
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION3D_H_ */
