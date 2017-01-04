
#ifndef SRC_CONFIG_LINEARELASTICITY2D_H_
#define SRC_CONFIG_LINEARELASTICITY2D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

struct LinearElasticity2DMaterial: public Configuration {

	enum MODEL {
		LINEAR_ELASTIC_ISOTROPIC = 0,
		LINEAR_ELASTIC_ORTHOTROPIC = 1,
		LINEAR_ELASTIC_ANISOTROPIC = 2
	};

	enum Parameter {
		DENSITY = 0,
		HEAT_CAPACITY,
		POISSON_RATIO,
		YOUNG_MODULUS_X,
		YOUNG_MODULUS_Y,
		THERMAL_EXPANSION_X,
		THERMAL_EXPANSION_Y
	};

	PARAMETER(std::string, density, "Density"                , "0");
	PARAMETER(std::string, Cp     , "Termal capacity."       , "0");
	PARAMETER(std::string, MI     , "Poisson ratio."         , "0.3");
	PARAMETER(std::string, EX     , "Young modulus X."       , "1");
	PARAMETER(std::string, EY     , "Young modulus Y."       , "1");
	PARAMETER(std::string, TEX    , "Thermal expansion X."   , "1");
	PARAMETER(std::string, TEY    , "Thermal expansion Y."   , "1");

	OPTION(MODEL, model, "Material model", MODEL::LINEAR_ELASTIC_ISOTROPIC, OPTIONS({
		{ "LINEAR_ELASTIC_ISOTROPIC"  , MODEL::LINEAR_ELASTIC_ISOTROPIC  , "Isotropic material." },
		{ "LINEAR_ELASTIC_ORTHOTROPIC", MODEL::LINEAR_ELASTIC_ORTHOTROPIC, "Orthotropic material." },
		{ "LINEAR_ELASTIC_ANISOTROPIC", MODEL::LINEAR_ELASTIC_ANISOTROPIC, "Anisotropic material." }
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
