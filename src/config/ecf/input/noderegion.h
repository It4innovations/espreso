
#ifndef SRC_CONFIG_ECF_INPUT_NODEREGION_H_
#define SRC_CONFIG_ECF_INPUT_NODEREGION_H_

#include "config/description.h"

#include <string>

namespace espreso {

struct InputNodeRegionConfiguration: public ECFDescription {

	enum class SHAPE {
		CIRCLE,
		BLOCK
	};

	SHAPE shape;
	double x, y, z, lx, ly, lz, radius;

	InputNodeRegionConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_INPUT_NODEREGION_H_ */
