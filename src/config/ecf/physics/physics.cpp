
#include "physics.h"
#include "config/configuration.hpp"
#include "config/conditions.h"

espreso::PhysicsConfiguration::PhysicsConfiguration(DIMENSION dim, MaterialConfiguration::PHYSICAL_MODEL physicalModel)
: dimension(dim), physical_model(physicalModel)
{
	load_steps = 1;
	REGISTER(load_steps, ECFMetaData()
			.setdescription({ "Number of loadSteps" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	ecfdescription->addSpace();

	interpolation = INTERPOLATION::LINEAR;
	REGISTER(interpolation, ECFMetaData()
			.setdescription({ "Data interpolation" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear interpolation."))
			.addoption(ECFOption().setname("QUADRATIC").setdescription("Quadratic interpolation.")));

	REGISTER(discretization, ECFMetaData()
		.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::OPTION })
		.setdescription({ "Discretization settings for regions.", "Discretization of stiffness matrices" })
		.setpattern({ "MY_REGION", "FEM_LOADED" })
		.addoption(ECFOption().setname("FEM_LOADED").setdescription("Finite elements from loaded mesh."))
		.addoption(ECFOption().setname("FEM_LINEAR").setdescription("Finite elements always linear."))
		.addoption(ECFOption().setname("FEM_QUADRATIC").setdescription("Finite elements always quadratic."))
		.addoption(ECFOption().setname("FEM_TDNNS").setdescription("Finite elements TDNNS."))
		.addoption(ECFOption().setname("BEM").setdescription("Boundary elements.")));

	ecfdescription->addSeparator();

	REGISTER(materials, ECFMetaData()
			.setdescription({ "The name of a material.", "Material description" })
			.setdatatype({ ECFDataType::STRING })
			.setpattern({ "MY_MATERIAL" }),
			&dimension, physicalModel);

	REGISTER(material_set, ECFMetaData()
			.setdescription({ "The name of a region.", "The name of a material." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::MATERIAL })
			.setpattern({ "MY_REGION", "MY_MATERIAL" }));

	ecfdescription->addSeparator();

	contact_interfaces = false;
	REGISTER(contact_interfaces, ECFMetaData()
			.setdescription({ "Consistent stabilization" })
			.setdatatype({ ECFDataType::BOOL }));

	ecfdescription->addSeparator();

	REGISTER(initial_temperature, ECFMetaData()
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setdescription({ "The name of a region.", "Initial temperature" })
			.setpattern({ "MY_REGION", "273.15" }),
			ECFMetaData().getboundaryconditionvariables(), "273.15");

	if (dimension == DIMENSION::D2) {
		REGISTER(thickness, ECFMetaData()
				.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
				.setdescription({ "The name of a region.", "Thickness" })
				.setpattern({ "MY_REGION", "1" }),
				ECFMetaData().getboundaryconditionvariables(), "1");
	}
}
