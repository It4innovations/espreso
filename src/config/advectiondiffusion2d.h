

#ifndef SRC_CONFIG_ADVECTIONDIFFUSION2D_H_
#define SRC_CONFIG_ADVECTIONDIFFUSION2D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

struct AdvectionDiffusion2DConfiguration: public Configuration {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	enum class HEAT_FLUX {
		GENERAL = 0,
		CONVECTIVE = 1
	};

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");

	OPTION(STABILIZATION, stabilization, "The type of the stabilization.", STABILIZATION::SUPG, OPTIONS({
		{ "SUPG", STABILIZATION::SUPG, "SUPG stabilization." },
		{ "CAU" , STABILIZATION::CAU , "CAU stabilization." },
	}));
	PARAMETER(double, sigma, "Inconsistent stabilization parameters.", 0);

	OPTION(HEAT_FLUX, heat_flux, "The type of heat flux condition.", HEAT_FLUX::GENERAL, OPTIONS({
		{ "GENERAL"    , HEAT_FLUX::GENERAL   , "General heat flux." },
		{ "CONVECTIVE" , HEAT_FLUX::CONVECTIVE, "Convective heat flux with parameters heat transfer and external temperature." },
	}));

	SUBMAP(std::string, std::string, general_heat_flux        , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, heat_transfer_coefficient, "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, external_temperature     , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBMAP(std::string, std::string, initial_temperature , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, temperature         , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, heat_source         , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, translation_motions , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, thickness           , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTOR(MaterialParameters, materials, "Vector of materials (counterd from 1).", "1", "Description of material with index 1");
	SUBMAP(size_t, std::string  , material_set, "Assign materials to regions", "<MATERIAL_INDEX>", "<REGION>");
};

}



#endif /* SRC_CONFIG_ADVECTIONDIFFUSION2D_H_ */
