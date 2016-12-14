
#ifndef SRC_CONFIG_RESULTS_H_
#define SRC_CONFIG_RESULTS_H_

#include "configuration.h"

namespace espreso {

struct Results: public Configuration {

	PARAMETER(double, norm, "Norm of the solution", 0);
};

}




#endif /* SRC_CONFIG_RESULTS_H_ */
