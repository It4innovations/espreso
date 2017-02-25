
#ifndef SRC_CONFIGURATION_MATERIAL_COORDINATESYSTEM_H_
#define SRC_CONFIGURATION_MATERIAL_COORDINATESYSTEM_H_

#include "../configuration.hpp"

namespace espreso {

struct CoordinateSystem: public Configuration {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	OPTION(TYPE, type, "System type", TYPE::CARTESIAN, OPTIONS({
		{ "CARTESIAN"  , TYPE::CARTESIAN  , "Cartesian system [CENTER, ROTATION]." },
		{ "CYLINDRICAL", TYPE::CYLINDRICAL, "Cylindrical system [CENTER, ROTATION]." },
		{ "SPHERICAL"  , TYPE::SPHERICAL  , "Spherical system [CENTER]." }
	}));

	PARAMETER(std::string, rotation_x, "Rotation 'x' of Cartesian system (in angles).", "0");
	PARAMETER(std::string, rotation_y, "Rotation 'y' of Cartesian system (in angles).", "0");
	PARAMETER(std::string, rotation_z, "Rotation 'z' of Cartesian system (in angles).", "0");

	PARAMETER(std::string, center_x, "x-center.", "0");
	PARAMETER(std::string, center_y, "x-center.", "0");
	PARAMETER(std::string, center_z, "x-center.", "0");
};


}



#endif /* SRC_CONFIGURATION_MATERIAL_COORDINATESYSTEM_H_ */
