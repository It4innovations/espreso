
#ifndef SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS3D_H_
#define SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS3D_H_

#include "structuralmechanics.h"
#include "discretization.h"

namespace espreso {

struct StructuralMechanics3DMaterial: public Configuration {

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::DENSITY>            , density, "Density"             , {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::HEAT_CAPACITY>      , Cp     , "Termal capacity."    , {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::POISSON_RATIO_XY>   , MIXY   , "Poisoon ratio XY."   , {"0.3"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::POISSON_RATIO_XZ>   , MIXZ   , "Poisoon ratio XZ."   , {"0.3"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::POISSON_RATIO_YZ>   , MIYZ   , "Poisoon ratio YZ."   , {"0.3"});

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::YOUNG_MODULUS_X>    , EX     , "Young modulus X."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::YOUNG_MODULUS_Y>    , EY     , "Young modulus Y."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::YOUNG_MODULUS_Z>    , EZ     , "Young modulus Z."    , {"1"});

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_EXPANSION_X>, TEX    , "Thermal expansion X.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_EXPANSION_Y>, TEY    , "Thermal expansion Y.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_EXPANSION_Z>, TEZ    , "Thermal expansion Z.", {"1"});

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::SHEAR_MODULUS_XY>   , GXY    , "Shear modulus XY."   , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::SHEAR_MODULUS_XZ>   , GXZ    , "Shear modulus XZ."   , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::SHEAR_MODULUS_YZ>   , GYZ    , "Shear modulus YZ."   , {"1"});

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D11>                , D11    , "Coefficient D11."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D12>                , D12    , "Coefficient D12."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D13>                , D13    , "Coefficient D13."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D14>                , D14    , "Coefficient D14."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D15>                , D15    , "Coefficient D15."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D16>                , D16    , "Coefficient D16."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D22>                , D22    , "Coefficient D22."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D23>                , D23    , "Coefficient D23."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D24>                , D24    , "Coefficient D24."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D25>                , D25    , "Coefficient D25."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D26>                , D26    , "Coefficient D26."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D33>                , D33    , "Coefficient D33."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D34>                , D34    , "Coefficient D34."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D35>                , D35    , "Coefficient D35."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D36>                , D36    , "Coefficient D36."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D44>                , D44    , "Coefficient D44."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D45>                , D45    , "Coefficient D45."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D46>                , D46    , "Coefficient D46."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D55>                , D55    , "Coefficient D55."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D56>                , D56    , "Coefficient D56."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::D66>                , D66    , "Coefficient D66."    , {"1"});

	OPTION(MATERIAL_MODEL, model, "Material model", MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC, OPTIONS({
		{ "LINEAR_ELASTIC_ISOTROPIC"  , MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC  , "Isotropic material." },
		{ "LINEAR_ELASTIC_ORTHOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC, "Orthotropic material." },
		{ "LINEAR_ELASTIC_ANISOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC, "Anisotropic material." }
	}));

	SUBCONFIG(CoordinateSystem, coordinate_system, "Element coordinate system.");
};

struct StructuralMechanics3DConfiguration: public StructuralMechanicsConfiguration {

	SUBMAPTOCONFIG(std::string, StructuralMechanics3DMaterial, materials, "Material description.", "MATERIAL", "Material description");

	OPTION(DISCRETIZATION, discretization, "Discretization type", DISCRETIZATION::FEM, OPTIONS({
		{ "FEM", DISCRETIZATION::FEM, "Finite element method." },
		{ "BEM", DISCRETIZATION::BEM, "Boundary element method." }
	}));
};

}


#endif /* SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS3D_H_ */
