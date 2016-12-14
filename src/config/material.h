
#ifndef SRC_CONFIG_MATERIAL_H_
#define SRC_CONFIG_MATERIAL_H_

#include "configuration.h"

namespace espreso {

enum class MATERIAL_MODEL {
	LINEAR_ELASTIC_ISOTROPIC = 0,
	LINEAR_ELASTIC_ORTHOTROPIC = 1,
	LINEAR_ELASTIC_ANISOTROPIC = 2
};

struct MaterialParameters: public Configuration {

	PARAMETER(std::string, DENS, "Density"             , "7850");
	PARAMETER(std::string, MIXY, "Poisson ratio XY."   , "0.3");
	PARAMETER(std::string, MIXZ, "Poisson ratio XZ."   , "0.3");
	PARAMETER(std::string, MIYZ, "Poisson ratio YZ."   , "0.3");
	PARAMETER(std::string, EX  , "Young modulus X."    , "2.1e11");
	PARAMETER(std::string, EY  , "Young modulus Y."    , "2.1e11");
	PARAMETER(std::string, EZ  , "Young modulus Z."    , "2.1e11");
	PARAMETER(std::string, C   , "Termal capacity."    , "1");
	PARAMETER(std::string, KXX , "Termal conduction X.", "1");
	PARAMETER(std::string, KXY , "Termal conduction Y.", "1");
	PARAMETER(std::string, KXZ , "Termal conduction Z.", "1");
	PARAMETER(std::string, ALPX, "Termal expansion X." , "1");
	PARAMETER(std::string, ALPY, "Termal expansion Y." , "1");
	PARAMETER(std::string, ALPZ, "Termal expansion Z." , "1");

	OPTION(MATERIAL_MODEL, model, "Material model", MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC, OPTIONS({
		{ "LINEAR_ELASTIC_ISOTROPIC"  , MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC  , "Isotropic." },
		{ "LINEAR_ELASTIC_ORTHOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC, "Orthotropic." },
		{ "LINEAR_ELASTIC_ANISOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC, "Anisotropic." }
	}));
};

}




#endif /* SRC_CONFIG_MATERIAL_H_ */
