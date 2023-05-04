
#include "plasticityproperties.h"
#include "config/configuration.hpp"

espreso::PlasticityPropertiesConfiguration::PlasticityPropertiesConfiguration(DIMENSION *D)
: dimension(D)
{
	REGISTER(initial_yield_stress, ECFMetaData()
			.setdescription({ "Initial yield stress." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(isotropic_hardening, ECFMetaData()
			.setdescription({ "Isotropic hardening." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(kinematic_hardening, ECFMetaData()
			.setdescription({ "Kinematic hardening." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}




