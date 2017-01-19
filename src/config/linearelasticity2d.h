
#ifndef SRC_CONFIG_LINEARELASTICITY2D_H_
#define SRC_CONFIG_LINEARELASTICITY2D_H_

#include "solver.h"
#include "materialholder.h"

namespace espreso {

struct LinearElasticity2DMaterial: public Configuration {

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
};

struct LinearElasticity2DConfiguration: public Configuration {

	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");

	OPTION(ELEMENT_BEHAVIOUR, element_behaviour, "The type elements.", ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS, OPTIONS({
		{ "PLAIN_STRAIN"               , ELEMENT_BEHAVIOUR::PLANE_STRAIN, "Strain element." },
		{ "AXISYMMETRIC"               , ELEMENT_BEHAVIOUR::AXISYMMETRIC, "Axisymmetric element." },
		{ "PLANE_STRESS"               , ELEMENT_BEHAVIOUR::PLANE_STRESS, "Stress element." },
		{ "PLANE_STRESS_WITH_THICKNESS", ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS , "Stress element with thickness." },
	}));
	PARAMETER(double, angular_velocity_x, "Angular velocity for x-axis.", 0);
	PARAMETER(double, angular_velocity_y, "Angular velocity for x-axis.", 0);
	PARAMETER(double, angular_velocity_z, "Angular velocity for x-axis.", 0);

	SUBMAP(std::string, std::string, displacement       , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, normal_presure     , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, initial_temperature, "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, temperature        , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, acceleration       , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, obstacle           , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, normal_direction   , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, thickness          , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTOR(LinearElasticity2DMaterial, materials, "Vector of materials.", "1", "Description of material '1'");
	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "<REGION>", "<MATERIAL_NAME>");

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);
};

}




#endif /* SRC_CONFIG_LINEARELASTICITY2D_H_ */
