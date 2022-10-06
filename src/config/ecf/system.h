
#ifndef SRC_CONFIG_ECF_SYSTEM_H_
#define SRC_CONFIG_ECF_SYSTEM_H_

#include "config/description.h"

namespace espreso {

struct SystemConfiguration: public ECFDescription {

	int mpi_affinity;

	SystemConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_SYSTEM_H_ */
