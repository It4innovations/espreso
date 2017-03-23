
#ifndef SRC_CONFIGURATION_INPUT_INPUTGENERATORGRIDTOWER_H_
#define SRC_CONFIGURATION_INPUT_INPUTGENERATORGRIDTOWER_H_

#include "inputgeneratorgrid.h"

namespace espreso {

struct NamedGridConfiguration: public GridConfiguration {

	PARAMETER(std::string, name, "Name of a body", "BODY");
};

struct GridTowerConfiguration: public Configuration {

	enum class DIRECTION {
		X,
		Y,
		Z
	};

	OPTION(DIRECTION, direction, "Direction of the tower", DIRECTION::X, OPTIONS({
		{ "X", DIRECTION::X, "Grids are placed at X direction." },
		{ "Y", DIRECTION::Y, "Grids are placed at Y direction." },
		{ "Z", DIRECTION::Z, "Grids are placed at Z direction." }
	}));

	SUBMAPTOCONFIG(size_t, NamedGridConfiguration, grids, "list of grids", "0", "Configuration of grid with index '0'");
};

}



#endif /* SRC_CONFIGURATION_INPUT_INPUTGENERATORGRIDTOWER_H_ */
