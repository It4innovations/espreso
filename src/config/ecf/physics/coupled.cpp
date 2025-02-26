
#include "coupled.h"
#include "config/configuration.hpp"

using namespace espreso;

ThermoElasticityLoadStepConfiguration::ThermoElasticityLoadStepConfiguration()
{
	coupling = COUPLING::WEAK;
	REGISTER(coupling, ECFMetaData()
			.setdescription({ "Coupling type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("WEAK").setdescription("Physics are computed separately."))
			.addoption(ECFOption().setname("STRONG").setdescription("Physics are computed at once")));

	REGISTER(heat_transfer, ECFMetaData()
			.setdescription({ "Heat transfer configuration" }));

	REGISTER(structural_mechanics, ECFMetaData()
			.setdescription({ "Structural mechanics configuration" }));
}

ThermoElasticityConfiguration::ThermoElasticityConfiguration()
: PhysicsConfiguration((MaterialConfiguration::PHYSICAL_MODEL)(MaterialConfiguration::PHYSICAL_MODEL::THERMAL | MaterialConfiguration::PHYSICAL_MODEL::STRUCTURAL_MECHANICS)),
  HeatTransferGlobalSettings(ecfdescription),
  StructuralMechanicsGlobalSettings(ecfdescription)
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Settings for each load step", "LoadStep" })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }));
}
