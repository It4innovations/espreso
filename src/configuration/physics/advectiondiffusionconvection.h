
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONCONVECTION_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONCONVECTION_H_

#include "../configuration.hpp"

namespace espreso {

struct AdvectionDiffusionConvection: public Configuration {

	PARAMETER(std::string, external_temperature     , "External temperature."     , "0");
	PARAMETER(std::string, heat_transfer_coefficient, "Heat transfer coefficient.", "0");
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONCONVECTION_H_ */
