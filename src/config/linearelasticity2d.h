
#ifndef SRC_CONFIG_LINEARELASTICITY2D_H_
#define SRC_CONFIG_LINEARELASTICITY2D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

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

	SUBVECTOR(MaterialParameters, materials   , "Vector of materials (counterd from 1).", "1", "Description of material with index 1");
	SUBMAP(size_t, std::string  , material_set, "Assign materials to regions", "<MATERIAL_INDEX>", "<REGION>");
};

}




#endif /* SRC_CONFIG_LINEARELASTICITY2D_H_ */
