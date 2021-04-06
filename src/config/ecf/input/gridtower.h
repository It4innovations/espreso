
#ifndef SRC_CONFIG_ECF_INPUT_GRIDTOWER_H_
#define SRC_CONFIG_ECF_INPUT_GRIDTOWER_H_

#include "grid.h"

namespace espreso {

struct GridTowerGeneratorConfiguration: public ECFDescription {

	enum class DIRECTION {
		X,
		Y,
		Z
	};

	DIRECTION direction;

	std::map<esint, GridGeneratorConfiguration> grids;

	GridTowerGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_GRIDTOWER_H_ */
