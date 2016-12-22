

#ifndef SRC_CONFIG_ADVECTIONDIFFUSION2D_H_
#define SRC_CONFIG_ADVECTIONDIFFUSION2D_H_

#include "solver.h"
#include "advectiondiffusionconvection.h"

namespace espreso {

struct AdvectionDiffusion2DMaterial: public Configuration {

	enum MODEL {
		ISOTROPIC = 0,
		DIAGONAL = 1,
		SYMMETRIC = 2,
		ANISOTROPIC = 3
	};

	enum Parameter {
		DENSITY = 0,
		HEAT_CAPACITY,
		THERMAL_CONDUCTIVITY_XX,
		THERMAL_CONDUCTIVITY_YY,
		THERMAL_CONDUCTIVITY_XY,
		THERMAL_CONDUCTIVITY_YX
	};

	PARAMETER(std::string, density, "Density"                , "0");
	PARAMETER(std::string, Cp     , "Termal capacity."       , "0");
	PARAMETER(std::string, KXX    , "Termal conductivity XX.", "1");
	PARAMETER(std::string, KYY    , "Termal conductivity YY.", "1");
	PARAMETER(std::string, KXY    , "Termal conductivity XY.", "1");
	PARAMETER(std::string, KYX    , "Termal conductivity YX.", "1");

	OPTION(MODEL, model, "Material model", MODEL::ISOTROPIC, OPTIONS({
		{ "ISOTROPIC"  , MODEL::ISOTROPIC  , "Isotropic." },
		{ "DIAGONAL"   , MODEL::DIAGONAL   , "Diagonal." },
		{ "SYMMETRIC"  , MODEL::SYMMETRIC  , "Symmetric." },
		{ "ANISOTROPIC", MODEL::ANISOTROPIC, "Anisotropic." }
	}));
};

struct AdvectionDiffusion2DConfiguration: public Configuration {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
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

	SUBMAP(std::string, std::string, heat_flux                , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, heat_flow                , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTOR(AdvectionDiffusionConvection, convection, "Region with convective heat flux", "<REGION>", "Convective parameters.");

	SUBMAP(std::string, std::string, initial_temperature , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, temperature         , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, heat_source         , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, translation_motions , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, thickness           , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTOR(AdvectionDiffusion2DMaterial, materials, "Vector of materials.", "1", "Description of material '1'");
	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "<MATERIAL_NAME>", "<REGION>");

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);
};

}



#endif /* SRC_CONFIG_ADVECTIONDIFFUSION2D_H_ */
