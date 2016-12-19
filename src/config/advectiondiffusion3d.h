
#ifndef SRC_CONFIG_ADVECTIONDIFFUSION3D_H_
#define SRC_CONFIG_ADVECTIONDIFFUSION3D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

struct AdvectionDiffusion3DMaterial: public Configuration {

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
		THERMAL_CONDUCTIVITY_ZZ,
		THERMAL_CONDUCTIVITY_XY,
		THERMAL_CONDUCTIVITY_XZ,
		THERMAL_CONDUCTIVITY_YZ,
		THERMAL_CONDUCTIVITY_YX,
		THERMAL_CONDUCTIVITY_ZX,
		THERMAL_CONDUCTIVITY_ZY
	};

	PARAMETER(std::string, density, "Density"                , "0");
	PARAMETER(std::string, Cp     , "Termal capacity."       , "0");
	PARAMETER(std::string, KXX    , "Termal conductivity XX.", "1");
	PARAMETER(std::string, KYY    , "Termal conductivity YY.", "1");
	PARAMETER(std::string, KZZ    , "Termal conductivity ZZ.", "1");

	PARAMETER(std::string, KXY    , "Termal conductivity XY.", "1");
	PARAMETER(std::string, KXZ    , "Termal conductivity XZ.", "1");
	PARAMETER(std::string, KYZ    , "Termal conductivity YZ.", "1");

	PARAMETER(std::string, KYX    , "Termal conductivity YX.", "1");
	PARAMETER(std::string, KZX    , "Termal conductivity ZX.", "1");
	PARAMETER(std::string, KZY    , "Termal conductivity ZY.", "1");

	OPTION(MODEL, model, "Material model", MODEL::ISOTROPIC, OPTIONS({
		{ "ISOTROPIC"  , MODEL::ISOTROPIC  , "Isotropic." },
		{ "DIAGONAL"   , MODEL::DIAGONAL   , "Diagonal." },
		{ "SYMMETRIC"  , MODEL::SYMMETRIC  , "Symmetric." },
		{ "ANISOTROPIC", MODEL::ANISOTROPIC, "Anisotropic." }
	}));
};

struct AdvectionDiffusion3DConfiguration: public Configuration {

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");
};

}



#endif /* SRC_CONFIG_ADVECTIONDIFFUSION3D_H_ */
