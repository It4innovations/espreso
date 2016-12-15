
#ifndef SRC_CONFIG_ADVECTIONDIFFUSION3D_H_
#define SRC_CONFIG_ADVECTIONDIFFUSION3D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

enum class MATERIAL_MODEL_AD3D {
	ISOTROPIC = 0,
	DIAGONAL = 1,
	SYMMETRIC = 2,
	ANISOTROPIC = 3
};

struct AdvectionDiffusion3DMaterial: public Configuration {

	PARAMETER(std::string, DENS, "Density"                , "7850");
	PARAMETER(std::string, Cp  , "Termal capacity."       , "1");
	PARAMETER(std::string, KXX , "Termal conductivity XX.", "1");
	PARAMETER(std::string, KYY , "Termal conductivity YY.", "1");
	PARAMETER(std::string, KZZ , "Termal conductivity YY.", "1");
	PARAMETER(std::string, KXY , "Termal conductivity XY.", "1");
	PARAMETER(std::string, KXZ , "Termal conductivity XZ.", "1");
	PARAMETER(std::string, KYX , "Termal conductivity YX.", "1");
	PARAMETER(std::string, KYZ , "Termal conductivity YZ.", "1");
	PARAMETER(std::string, KZX , "Termal conductivity ZX.", "1");
	PARAMETER(std::string, KZY , "Termal conductivity ZY.", "1");

	OPTION(MATERIAL_MODEL_AD3D, model, "Material model", MATERIAL_MODEL_AD3D::ISOTROPIC, OPTIONS({
		{ "ISOTROPIC"  , MATERIAL_MODEL_AD3D::ISOTROPIC  , "Isotropic." },
		{ "DIAGONAL"   , MATERIAL_MODEL_AD3D::DIAGONAL   , "Diagonal." },
		{ "SYMMETRIC"  , MATERIAL_MODEL_AD3D::SYMMETRIC  , "Symmetric." },
		{ "ANISOTROPIC", MATERIAL_MODEL_AD3D::ANISOTROPIC, "Anisotropic." }
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
