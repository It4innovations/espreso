
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONRADIATION_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONRADIATION_H_

#include "../configuration.hpp"

namespace espreso {


struct AdvectionDiffusionRadiation: public Configuration {

	PARAMETER(std::string, emissivity          , "Surface emissivity."  , "0");
	PARAMETER(std::string, external_temperature, "External temperature.", "0");
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONRADIATION_H_ */
