
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
{
	REGISTER(impedance, ECFMetaData()
			.setdescription({ "Impedance value" })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

PointSourceConfiguration::PointSourceConfiguration()
{
	type = TYPE::USER;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Point source type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("USER").setdescription("User defined."))
			.addoption(ECFOption().setname("POWER").setdescription("Power."))
			.addoption(ECFOption().setname("INTENSITY").setdescription("Intensity."))
			.addoption(ECFOption().setname("FLOW").setdescription("Flow.")));

	REGISTER(volume_flow, ECFMetaData()
			.setdescription({ "Volume flow rate" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(phase, ECFMetaData()
			.setdescription({ "Phase" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(reference_intensity, ECFMetaData()
			.setdescription({ "Free space reference intensity" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(distance_from_source, ECFMetaData()
			.setdescription({ "Distance from center" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(reference_power, ECFMetaData()
			.setdescription({ "Free space reference power" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(source_amplitude, ECFMetaData()
			.setdescription({ "Source amplitude" })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

AcousticLoadStepConfiguration::AcousticLoadStepConfiguration(DIMENSION *D)
{
	system = SYSTEM::REAL;
	REGISTER(system, ECFMetaData()
			.setdescription({ "Linear system type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("REAL").setdescription("Solve with real values."))
			.addoption(ECFOption().setname("COMPLEX").setdescription("Solve with complex values.")));

	REGISTER(monopole_source, ECFMetaData()
			.setdescription({ "The name of a region.", "Monopole source" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" })
			.setdynamic());

	REGISTER(dipole_source, ECFMetaData()
			.setdescription({ "The name of a region.", "Dipole source." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" }),
			D);

	REGISTER(acoustic_pressure, ECFMetaData()
			.setdescription({ "The name of a region.", "Pressure" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" })
			.setdynamic());

	REGISTER(normal_acceleration, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal Acceleration" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" })
			.setdynamic());

	REGISTER(acceleration, ECFMetaData()
			.setdescription({ "The name of a region.", "Acceleration" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" })
			.setdynamic(),
			D);

	REGISTER(impedance, ECFMetaData()
			.setdescription({ "The name of a region.", "Impedance" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION"})
			.setdynamic());

	REGISTER(point_source, ECFMetaData()
			.setdescription({ "The name of a region.", "Point source" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" })
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
