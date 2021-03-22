
#include <config/ecf/ecf.h>
#include "config/configuration.hpp"

using namespace espreso;

AcousticsLoadStepConfiguration::AcousticsLoadStepConfiguration(DIMENSION *D)
{
	REGISTER(acoustic_pressure, ECFMetaData()
			.setdescription({ "The name of a region.", "Pressure" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" })
			.setdynamic(),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(normal_acceleration, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal Acceleration" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" })
			.setdynamic(),
			ECFMetaData::getboundaryconditionvariables());
}

AcousticsConfiguration::AcousticsConfiguration(DIMENSION d)
: PhysicsConfiguration(d, MaterialConfiguration::PHYSICAL_MODEL::ACOUSTICS),
  dimension(d)
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Settings for each load step.", "Load step index." })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }),
			&dimension);
}