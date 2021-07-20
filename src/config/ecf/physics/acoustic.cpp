
#include <config/ecf/ecf.h>
#include "config/configuration.hpp"

using namespace espreso;

AcousticGlobalSettings::AcousticGlobalSettings(ECFObject *ecfdescription)
{
	init_temp_respect_bc = true;
	REGISTER(init_temp_respect_bc, ECFMetaData()
			.setdescription({ "Initial temperature follows BC" })
			.setdatatype({ ECFDataType::BOOL })
			.setform());

	diffusion_split = false;
	REGISTER(diffusion_split, ECFMetaData()
			.setdescription({ "Thermal shock stabilization" })
			.setdatatype({ ECFDataType::BOOL })
			.setform());
}

AcousticLoadStepConfiguration::AcousticLoadStepConfiguration(DIMENSION *D)
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

AcousticConfiguration::AcousticConfiguration(DIMENSION d)
: PhysicsConfiguration(d, MaterialConfiguration::PHYSICAL_MODEL::ACOUSTIC),
  AcousticGlobalSettings(ecfdescription),
  dimension(d)
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Settings for each load step.", "Load step index." })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }),
			&dimension);
}
