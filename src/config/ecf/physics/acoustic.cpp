
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

ImpedanceConfiguration::ImpedanceConfiguration()
: impedance(ECFExpression::Scope::BGPS)
{
	REGISTER(impedance, ECFMetaData()
			.setdescription({ "Impedance value" })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

AcousticLoadStepConfiguration::AcousticLoadStepConfiguration(DIMENSION *D)
{
	REGISTER(acoustic_pressure, ECFMetaData()
			.setdescription({ "The name of a region.", "Pressure" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" })
			.setdynamic(),
			ECFExpression::Scope::NODE);
	REGISTER(normal_acceleration, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal Acceleration" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" })
			.setdynamic(),
			ECFExpression::Scope::BGPS);

	REGISTER(impedance, ECFMetaData()
			.setdescription({ "The name of a region.", "Impedance" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION"})
			.setdynamic());
}

AcousticConfiguration::AcousticConfiguration(DIMENSION d)
: PhysicsConfiguration(d, MaterialConfiguration::PHYSICAL_MODEL::ACOUSTICS),
  AcousticGlobalSettings(ecfdescription),
  dimension(d)
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Settings for each load step.", "Load step index." })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }),
			&dimension);
}
