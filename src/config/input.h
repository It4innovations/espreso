
#ifndef SRC_CONFIG_INPUT_H_
#define SRC_CONFIG_INPUT_H_

#include "configuration.h"

namespace espreso {

struct ESPRESOInput: public Configuration {

	PARAMETER(std::string, path, "Path to an input description.", "");
	PARAMETER(size_t, domains, "Number of sub-domains of each cluster.", 8);
};

}



#endif /* SRC_CONFIG_INPUT_H_ */
