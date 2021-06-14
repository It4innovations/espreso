
#ifndef SRC_CONFIG_HOLDERS_RANGE_H_
#define SRC_CONFIG_HOLDERS_RANGE_H_

#include "config/configuration.h"

namespace espreso {

struct ECFRange {
	ECFParameter *parameter;

	std::string min, max, step;

	ECFRange();
};

}

#endif /* SRC_CONFIG_HOLDERS_RANGE_H_ */
