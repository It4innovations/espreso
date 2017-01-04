
#ifndef SRC_CONFIG_ADVECTIONDIFFUSION3D_H_
#define SRC_CONFIG_ADVECTIONDIFFUSION3D_H_

#include "solver.h"
#include "coordinatesystem.h"
#include "advectiondiffusionconvection.h"
#include "advectiondiffusionsolver.h"

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

	SUBCONFIG(CoordinateSystem, coordinate_system, "Element coordinate system.");
};

struct AdvectionDiffusion3DConfiguration: public Configuration {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	OPTION(STABILIZATION, stabilization, "The type of the stabilization.", STABILIZATION::SUPG, OPTIONS({
		{ "SUPG", STABILIZATION::SUPG, "SUPG stabilization." },
		{ "CAU" , STABILIZATION::CAU , "CAU stabilization." },
	}));
	PARAMETER(double, sigma, "Inconsistent stabilization parameters.", 0);

	SUBCONFIG(AdvectionDiffusionSolver, physics_solver, "Settings of physics solver.");

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");

	SUBMAP(std::string, std::string, heat_flux                , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, heat_flow                , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTOR(AdvectionDiffusionConvection, convection, "Region with convective heat flux", "<REGION>", "Convective parameters.");

	SUBMAP(std::string, std::string, initial_temperature , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, temperature         , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, heat_source         , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, translation_motions , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");

	SUBVECTOR(AdvectionDiffusion3DMaterial, materials, "Vector of materials.", "1", "Description of material '1'");
	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "<REGION>", "<MATERIAL_NAME>");

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);
};

}



#endif /* SRC_CONFIG_ADVECTIONDIFFUSION3D_H_ */
