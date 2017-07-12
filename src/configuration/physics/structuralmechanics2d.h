
#ifndef SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS2D_H_
#define SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS2D_H_

#include "structuralmechanics.h"

namespace espreso {

struct StructuralMechanics2DMaterial: public Configuration {

	PARAMETER(MaterialParam<MATERIAL_PARAMETER::DENSITY>            , density, "Density"             , {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::HEAT_CAPACITY>      , Cp     , "Termal capacity."    , {"0"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::POISSON_RATIO_XY>   , MI     , "Poisoon ratio."      , {"0.3"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::YOUNG_MODULUS_X>    , EX     , "Young modulus X."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::YOUNG_MODULUS_Y>    , EY     , "Young modulus Y."    , {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_EXPANSION_X>, TEX    , "Thermal expansion X.", {"1"});
	PARAMETER(MaterialParam<MATERIAL_PARAMETER::THERMAL_EXPANSION_Y>, TEY    , "Thermal expansion Y.", {"1"});

	OPTION(MATERIAL_MODEL, model, "Material model", MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC, OPTIONS({
		{ "LINEAR_ELASTIC_ISOTROPIC"  , MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC  , "Isotropic material." },
		{ "LINEAR_ELASTIC_ORTHOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC, "Orthotropic material." },
		{ "LINEAR_ELASTIC_ANISOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC, "Anisotropic material." }
	}));

	SUBCONFIG(CoordinateSystem, coordinate_system, "Element coordinate system.");
};

struct StructuralMechanics2DConfiguration: public StructuralMechanicsConfiguration {

	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	OPTION(ELEMENT_BEHAVIOUR, element_behaviour, "The type elements.", ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS, OPTIONS({
		{ "PLAIN_STRAIN"               , ELEMENT_BEHAVIOUR::PLANE_STRAIN, "Strain element." },
		{ "AXISYMMETRIC"               , ELEMENT_BEHAVIOUR::AXISYMMETRIC, "Axisymmetric element." },
		{ "PLANE_STRESS"               , ELEMENT_BEHAVIOUR::PLANE_STRESS, "Stress element." },
		{ "PLANE_STRESS_WITH_THICKNESS", ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS , "Stress element with thickness." },
	}));

	SUBMAPTOMAP(size_t, std::string, std::string, angular_velocity, "Angular velocity", "TIME_STEP", "Angular velocity settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, thickness, "Thickness", "TIME_STEP", "Thickness settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOCONFIG(std::string, StructuralMechanics2DMaterial, materials, "Material description.", "MATERIAL", "Material description");
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS2D_H_ */
