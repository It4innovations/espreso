
#include "coordinatesystem.h"

#include "config/conditions.h"
#include "config/configuration.hpp"

using namespace espreso;

CoordinateSystemConfiguration::CoordinateSystemConfiguration(DIMENSION *D)
: dimension(D), rotation(dimension, ECFMetaData::getcoordinatevariables(), "0"), center(dimension, ECFMetaData::getcoordinatevariables(), "0")
{
	type = TYPE::CARTESIAN;

	REGISTER(dimension, ECFMetaData()
				.setdescription({"Dimension"})
				.setdatatype({ ECFDataType::OPTION })
				.addoption(ECFOption().setname("D1").setdescription("D1"))
				.addoption(ECFOption().setname("D2").setdescription("D2"))
				.addoption(ECFOption().setname("D3").setdescription("D3"))
				.addoption(ECFOption().setname("Z").setdescription("Z"))
				.allowonly([&] () { return false; })
				.addconstraint(ECFFalseCondition()));

	REGISTER(type, ECFMetaData()
			.setdescription({ "Type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CARTESIAN").setdescription("Cartesian system."))
			.addoption(ECFOption().setname("CYLINDRICAL").setdescription("Cylindrical system."))
			.addoption(ECFOption().setname("SPHERICAL").setdescription("Spherical system.").allowonly([&] () { return *dimension == DIMENSION::D3; })));

	ecfdescription->addSeparator()->metadata.addconstraint(ECFCondition(type, ECFCondition::NOT_EQUALS, TYPE::CARTESIAN));

	REGISTER(center, ECFMetaData()
			.setname("Position")
			.setdescription({ "A center of the material." })
			.displayobjectname()
			.setgroup()
			.allowonly([&] () { return type != TYPE::CARTESIAN; })
			.addconstraint(ECFCondition(type, ECFCondition::NOT_EQUALS, TYPE::CARTESIAN)));
	
	ecfdescription->addSeparator()->metadata.addconstraint(ECFCondition(type, ECFCondition::EQUALS, TYPE::CARTESIAN));

	REGISTER(rotation, ECFMetaData()
			.setname("Rotation")
			.setdescription({ "A rotation of the material." })
			.displayobjectname()
			.setgroup()
			.allowonly([&] () { return type == TYPE::CARTESIAN; })
			.addconstraint(ECFCondition(type, ECFCondition::EQUALS, TYPE::CARTESIAN)));
}




