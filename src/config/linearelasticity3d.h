
#ifndef SRC_CONFIG_LINEARELASTICITY3D_H_
#define SRC_CONFIG_LINEARELASTICITY3D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

struct LinearElasticity3DConfiguration: public Configuration {

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");

	SUBMAP(std::string, std::string, displacement       , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, normal_presure     , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, initial_temperature, "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, temperature        , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, acceleration       , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, obstacle           , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, normal_direction   , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTOR(MaterialParameters, materials   , "Vector of materials (counterd from 1).", "1", "Description of material with index 1");
	SUBMAP(size_t, std::string  , material_set, "Assign materials to regions", "<MATERIAL_INDEX>", "<REGION>");
};

}



#endif /* SRC_CONFIG_LINEARELASTICITY3D_H_ */
