
#ifndef SRC_CONFIG_HOLDERS_RANGE_H_
#define SRC_CONFIG_HOLDERS_RANGE_H_

#include "config/configuration.h"

namespace espreso {

struct ECFRange {
	std::vector<ECFParameter*> parameter;

	std::string min, max, step;
};

}

#endif /* SRC_CONFIG_HOLDERS_RANGE_H_ */
