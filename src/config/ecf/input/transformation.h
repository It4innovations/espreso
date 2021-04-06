
#ifndef SRC_CONFIG_ECF_INPUT_TRANSFORMATION_H_
#define SRC_CONFIG_ECF_INPUT_TRANSFORMATION_H_

#include "config/description.h"

#include <string>

namespace espreso {

struct InputTransformationConfiguration: public ECFDescription {

	enum class TRANSFORMATION {
		TRANSLATE,
		ROTATE,
		SCALE,
		SHEAR
	};

	TRANSFORMATION transformation;
	double x, y, z;

	int instances;

	InputTransformationConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_TRANSFORMATION_H_ */
