

#ifndef SRC_CONFIGURATION_DECOMPOSER_H_
#define SRC_CONFIGURATION_DECOMPOSER_H_

#include "../configuration/configuration.hpp"

namespace espreso {

struct Decomposer: public Configuration {

	PARAMETER(std::string, parts, "Each MPI process will be decomposed into the specified number of parts (e.q. 1 2 4).", "1");
	PARAMETER(std::string, prefix, "Decomposition will be saved into PREFIX{PARTS} directories.", "DECOMPOSITION");
};

}



#endif /* SRC_CONFIGURATION_DECOMPOSER_H_ */
