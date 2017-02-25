
#ifndef SRC_CONFIGURATION_INPUT_INPUT_H_
#define SRC_CONFIGURATION_INPUT_INPUT_H_

#include "../configuration.hpp"

namespace espreso {

struct ESPRESOInput: public Configuration {

	PARAMETER(std::string, path, "Path to an input description.", "");
	PARAMETER(size_t, domains, "Number of sub-domains of each cluster.", 8);
};

}



#endif /* SRC_CONFIGURATION_INPUT_INPUT_H_ */
