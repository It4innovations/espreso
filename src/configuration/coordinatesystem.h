
#ifndef SRC_CONFIGURATION_COORDINATESYSTEM_H_
#define SRC_CONFIGURATION_COORDINATESYSTEM_H_

#include "../configuration/configuration.h"

namespace espreso {

struct CoordinateSystem: public Configuration {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	OPTION(TYPE, type, "System type", TYPE::CARTESIAN, OPTIONS({
		{ "CARTESIAN"  , TYPE::CARTESIAN  , "Cartesian system." },
		{ "CYLINDRICAL", TYPE::CYLINDRICAL, "Cylindrical system." },
		{ "SPHERICAL"  , TYPE::SPHERICAL  , "Spherical system." }
	}));

	PARAMETER(double, rotation_x, "Rotation 'x' of Cartesian system (in angles).", 0);
	PARAMETER(double, rotation_y, "Rotation 'y' of Cartesian system (in angles).", 0);
	PARAMETER(double, rotation_z, "Rotation 'z' of Cartesian system (in angles).", 0);

	PARAMETER(double, center_x, "x-center.", 0);
	PARAMETER(double, center_y, "x-center.", 0);
	PARAMETER(double, center_z, "x-center.", 0);
};


}



#endif /* SRC_CONFIGURATION_COORDINATESYSTEM_H_ */
