
#ifndef SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_
#define SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_

#include "basis/containers/point.h"

#include "config/holders/expression.h"
#include "config/description.h"

namespace espreso {

struct CoordinateSystemConfiguration: public ECFDescription {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	TYPE type;
	DIMENSION *dimension;

	ECFExpressionVector rotation;
	ECFExpressionVector center;

	CoordinateSystemConfiguration(DIMENSION *D);
};


}

#endif /* SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_ */
