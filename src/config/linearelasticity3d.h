
#ifndef SRC_CONFIG_LINEARELASTICITY3D_H_
#define SRC_CONFIG_LINEARELASTICITY3D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

struct LinearElasticity3DMaterial: public Configuration {

	enum MODEL {
		LINEAR_ELASTIC_ISOTROPIC = 0,
		LINEAR_ELASTIC_ORTHOTROPIC = 1,
		LINEAR_ELASTIC_ANISOTROPIC = 2
	};

	enum Parameter {
		DENSITY = 0,
		HEAT_CAPACITY,
		POISSON_RATIO_XY,
		POISSON_RATIO_XZ,
		POISSON_RATIO_YZ,
		YOUNG_MODULUS_X,
		YOUNG_MODULUS_Y,
		YOUNG_MODULUS_Z,
		THERMAL_EXPANSION_X,
		THERMAL_EXPANSION_Y,
		THERMAL_EXPANSION_Z,
		SHEAR_MODULUS_XY,
		SHEAR_MODULUS_XZ,
		SHEAR_MODULUS_YZ,

		D11,
		D12,
		D13,
		D14,
		D15,
		D16,
		D22,
		D23,
		D24,
		D25,
		D26,
		D33,
		D34,
		D35,
		D36,
		D44,
		D45,
		D46,
		D55,
		D56,
		D66
	};

	PARAMETER(std::string, density, "Density"             , "0");
	PARAMETER(std::string, Cp     , "Termal capacity."    , "0");
	PARAMETER(std::string, MIXY   , "Poisson ratio XY."   , "0.3");
	PARAMETER(std::string, MIXZ   , "Poisson ratio XZ."   , "0.3");
	PARAMETER(std::string, MIYZ   , "Poisson ratio YZ."   , "0.3");

	PARAMETER(std::string, EX     , "Young modulus X."    , "1");
	PARAMETER(std::string, EY     , "Young modulus Y."    , "1");
	PARAMETER(std::string, EZ     , "Young modulus Z."    , "1");

	PARAMETER(std::string, TEX    , "Thermal expansion X.", "1");
	PARAMETER(std::string, TEY    , "Thermal expansion Y.", "1");
	PARAMETER(std::string, TEZ    , "Thermal expansion Z.", "1");

	PARAMETER(std::string, GXY    , "Shear modulus XY."   , "1");
	PARAMETER(std::string, GXZ    , "Shear modulus XY."   , "1");
	PARAMETER(std::string, GYZ    , "Shear modulus YZ."   , "1");

	PARAMETER(std::string, d11    , "Coefficient."        , "1");
	PARAMETER(std::string, d12    , "Coefficient."        , "1");
	PARAMETER(std::string, d13    , "Coefficient."        , "1");
	PARAMETER(std::string, d14    , "Coefficient."        , "1");
	PARAMETER(std::string, d15    , "Coefficient."        , "1");
	PARAMETER(std::string, d16    , "Coefficient."        , "1");
	PARAMETER(std::string, d22    , "Coefficient."        , "1");
	PARAMETER(std::string, d23    , "Coefficient."        , "1");
	PARAMETER(std::string, d24    , "Coefficient."        , "1");
	PARAMETER(std::string, d25    , "Coefficient."        , "1");
	PARAMETER(std::string, d26    , "Coefficient."        , "1");
	PARAMETER(std::string, d33    , "Coefficient."        , "1");
	PARAMETER(std::string, d34    , "Coefficient."        , "1");
	PARAMETER(std::string, d35    , "Coefficient."        , "1");
	PARAMETER(std::string, d36    , "Coefficient."        , "1");
	PARAMETER(std::string, d44    , "Coefficient."        , "1");
	PARAMETER(std::string, d45    , "Coefficient."        , "1");
	PARAMETER(std::string, d46    , "Coefficient."        , "1");
	PARAMETER(std::string, d55    , "Coefficient."        , "1");
	PARAMETER(std::string, d56    , "Coefficient."        , "1");
	PARAMETER(std::string, d66    , "Coefficient."        , "1");

	OPTION(MODEL, model, "Material model", MODEL::LINEAR_ELASTIC_ISOTROPIC, OPTIONS({
		{ "LINEAR_ELASTIC_ISOTROPIC"  , MODEL::LINEAR_ELASTIC_ISOTROPIC  , "Isotropic material." },
		{ "LINEAR_ELASTIC_ORTHOTROPIC", MODEL::LINEAR_ELASTIC_ORTHOTROPIC, "Orthotropic material." },
		{ "LINEAR_ELASTIC_ANISOTROPIC", MODEL::LINEAR_ELASTIC_ANISOTROPIC, "Anisotropic material." }
	}));
};

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

	SUBVECTOR(LinearElasticity3DMaterial, materials   , "Vector of materials (counterd from 1).", "1", "Description of material with index 1");
	SUBMAP(size_t, std::string  , material_set, "Assign materials to regions", "<MATERIAL_INDEX>", "<REGION>");
};

}



#endif /* SRC_CONFIG_LINEARELASTICITY3D_H_ */
