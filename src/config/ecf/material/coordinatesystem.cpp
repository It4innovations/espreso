
#include "coordinatesystem.h"

#include "config/conditions.h"
#include "config/configuration.hpp"

using namespace espreso;

CoordinateSystemConfiguration::CoordinateSystemConfiguration(DIMENSION *D)
: dimension(D), rotation(dimension, ECFMetaData::getcoordinatevariables(), "0"), center(dimension, ECFMetaData::getcoordinatevariables(), "0")
{
	type = TYPE::CARTESIAN;

	REGISTER(type, ECFMetaData()
			.setdescription({ "Type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CARTESIAN").setdescription("Cartesian system."))
			.addoption(ECFOption().setname("CYLINDRICAL").setdescription("Cylindrical system."))
			.addoption(ECFOption().setname("SPHERICAL").setdescription("Spherical system.").allowonly([&] () { return *dimension == DIMENSION::D3; })));

	ecfdescription->addSeparator()->metadata.noexport();

	REGISTER(center, ECFMetaData()
			.setname("Position")
			.setdescription({ "A center of the material." })
			.displayobjectname()
			.setgroup()
			.allowonly([&] () { return type != TYPE::CARTESIAN; })
			.addconstraint(ECFCondition(type, ECFCondition::NOT_EQUALS, TYPE::CARTESIAN)));
	
	ecfdescription->addSeparator()->metadata.noexport();

	REGISTER(rotation, ECFMetaData()
			.setname("Rotation")
			.setdescription({ "A rotation of the material." })
			.displayobjectname()
			.setgroup()
			.allowonly([&] () { return type == TYPE::CARTESIAN; })
			.addconstraint(ECFCondition(type, ECFCondition::EQUALS, TYPE::CARTESIAN)));
}




