
#include "physics.h"
#include "config/configuration.hpp"
#include "config/conditions.h"

espreso::PhysicsConfiguration::PhysicsConfiguration(DIMENSION dim, MaterialConfiguration::PHYSICAL_MODEL physicalModel)
: dimension(dim), physical_model(physicalModel)
{
	load_steps = 1;
	REGISTER(load_steps, ECFMetaData()
			.setdescription({ "Number of loadSteps" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER })
			.setrange(1, 1000)
			.setform());

	ecfdescription->addSpace()->metadata.noexport();

	interpolation = INTERPOLATION::LINEAR;
	REGISTER(interpolation, ECFMetaData()
			.setdescription({ "Data interpolation" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear interpolation."))
			.addoption(ECFOption().setname("QUADRATIC").setdescription("Quadratic interpolation."))
			.setform());

	REGISTER(discretization, ECFMetaData()
		.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::OPTION })
		.setdescription({ "Discretization settings for regions.", "Discretization of stiffness matrices" })
		.setpattern({ "MY_REGION", "FEM_LOADED" })
		.setdynamic()
		.addoption(ECFOption().setname("FEM_LOADED").setdescription("Finite elements from loaded mesh."))
		.addoption(ECFOption().setname("FEM_LINEAR").setdescription("Finite elements always linear."))
		.addoption(ECFOption().setname("FEM_QUADRATIC").setdescription("Finite elements always quadratic."))
		.addoption(ECFOption().setname("FEM_TDNNS").setdescription("Finite elements TDNNS."))
		.addoption(ECFOption().setname("BEM").setdescription("Boundary elements.")));

	ecfdescription->addSeparator()->metadata.setform();

	REGISTER(materials, ECFMetaData()
			.setdescription({ "The name of a material.", "Material description" })
			.setdatatype({ ECFDataType::STRING })
			.setpattern({ "MY_MATERIAL" })
			.setdynamic()
			.setpatternname("Material"),
			&dimension, physicalModel);

	REGISTER(material_set, ECFMetaData()
			.setname("Material set")
			.setdescription({ "The name of a region.", "The name of a material." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::MATERIAL })
			.setpattern({ "MY_REGION", "MY_MATERIAL" })
			.setdynamic()
			.setpatternitemname("Material"));

	ecfdescription->addSeparator()->metadata.noexport();

	contact_interfaces = false;
	REGISTER(contact_interfaces, ECFMetaData()
			.setdescription({ "Consistent stabilization" })
			.setdatatype({ ECFDataType::BOOL })
			.setform());

	simd = true;
	REGISTER(simd, ECFMetaData()
			.setdescription({ "Assembler use intrinsic functions." })
			.setdatatype({ ECFDataType::BOOL })
			.setform());

	ecfdescription->addSeparator()->metadata.noexport();

	REGISTER(initial_temperature, ECFMetaData()
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setdescription({ "The name of a region.", "Initial temperature" })
			.setpattern({ "MY_REGION", "273.15" })
			.setdynamic(),
			"273.15", ECFExpression::Scope::EGPS);

	if (dimension == DIMENSION::D2) {
		REGISTER(thickness, ECFMetaData()
				.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
				.setdescription({ "The name of a region.", "Thickness" })
				.setpattern({ "MY_REGION", "1" })
				.setdynamic(),
				"1", ECFExpression::Scope::EGPS);
	}
}
