
#ifndef SRC_CONFIG_RESULTS_H_
#define SRC_CONFIG_RESULTS_H_

#include "configuration.h"

namespace espreso {

struct Results: public Configuration {

	PARAMETER(double, norm, "Norm of the solution", 0);

	PARAMETER(std::string, displacement_x, "x-displacement of each node.", "");
	PARAMETER(std::string, displacement_y, "y-displacement of each node.", "");
	PARAMETER(std::string, displacement_z, "z-displacement of each node.", "");
	PARAMETER(std::string, temperature   , "temperature of each node.", "");
	PARAMETER(std::string, pressure      , "presure of each node.", "");
};

}




#endif /* SRC_CONFIG_RESULTS_H_ */
