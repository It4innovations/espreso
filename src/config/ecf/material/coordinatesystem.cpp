
#include "coordinatesystem.h"

#include "config/conditions.h"
#include "config/configuration.hpp"

using namespace espreso;

CoordinateSystemConfiguration::CoordinateSystemConfiguration(DIMENSION *D)
: dimension(D), rotation(dimension), center(dimension)
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

bool CoordinateSystemConfiguration::isConst() const
{
	if (type == TYPE::CARTESIAN) {
		switch (*dimension) {
		case DIMENSION::D2: return rotation.z.evaluator->isConst();
		case DIMENSION::D3: return rotation.z.evaluator->isConst() && rotation.y.evaluator->isConst() && rotation.x.evaluator->isConst();
		}
	}
	switch (*dimension) {
	case DIMENSION::D2: return center.x.evaluator->isConst() && center.y.evaluator->isConst();
	case DIMENSION::D3: return center.x.evaluator->isConst() && center.y.evaluator->isConst() && center.z.evaluator->isConst();
	}
	return true;
}

bool CoordinateSystemConfiguration::isRotated() const
{
	if (type == TYPE::CARTESIAN) {
		switch (*dimension) {
		case DIMENSION::D2: return rotation.z.isset && rotation.z.evaluator->evaluate() != 0;
		case DIMENSION::D3:
			return
				(rotation.x.isset && rotation.x.evaluator->evaluate() != 0) |
				(rotation.y.isset && rotation.y.evaluator->evaluate() != 0) |
				(rotation.z.isset && rotation.z.evaluator->evaluate() != 0);
		default: return false;
		}
	}
	return true;
}



