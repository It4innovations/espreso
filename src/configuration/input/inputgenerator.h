
#ifndef SRC_CONFIGURATION_INPUT_INPUTGENERATOR_H_
#define SRC_CONFIGURATION_INPUT_INPUTGENERATOR_H_

#include "inputgeneratorgrid.h"
#include "inputgeneratorgridtower.h"
#include "inputgeneratorsphere.h"

namespace espreso {

enum class GENERATOR_SHAPE {
	GRID,
	GRID_TOWER,
	SPHERE
};

struct ESPRESOGenerator: public Configuration {

	OPTION(GENERATOR_SHAPE, shape, "Generated shape", GENERATOR_SHAPE::GRID, OPTIONS({
		{ "GRID"      , GENERATOR_SHAPE::GRID      , "Rectangular grid with empty spaces." },
		{ "GRID_TOWER", GENERATOR_SHAPE::GRID_TOWER, "Tower of rectangular grids." },
		{ "SPHERE"    , GENERATOR_SHAPE::SPHERE    , "Hollow sphere." }
	}));

	SUBCONFIG(GridConfiguration     , grid      , "Detailed specification of grid shape.");
	SUBCONFIG(GridTowerConfiguration, grid_tower, "Detailed specification of grid tower.");
	SUBCONFIG(SphereConfiguration   , sphere    , "Detailed specification of spherical shape.");
};

}



#endif /* SRC_CONFIGURATION_INPUT_INPUTGENERATOR_H_ */
