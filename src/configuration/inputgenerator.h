
#ifndef SRC_CONFIGURATION_INPUTGENERATOR_H_
#define SRC_CONFIGURATION_INPUTGENERATOR_H_

#include "../configuration/inputgeneratorgrid.h"
#include "../configuration/inputgeneratorsphere.h"

namespace espreso {

enum class GENERATOR_SHAPE {
	GRID,
	SPHERE
};

struct ESPRESOGenerator: public Configuration {

	OPTION(GENERATOR_SHAPE, shape, "Generated shape", GENERATOR_SHAPE::GRID, OPTIONS({
		{ "GRID"  , GENERATOR_SHAPE::GRID  , "Rectangular grid with empty spaces." },
		{ "SPHERE", GENERATOR_SHAPE::SPHERE, "Hollow sphere." }
	}));

	SUBCONFIG(GridConfiguration  , grid  , "Detailed specification of grid shape.");
	SUBCONFIG(SphereConfiguration, sphere, "Detailed specification of spherical shape.");
};

}



#endif /* SRC_CONFIGURATION_INPUTGENERATOR_H_ */
