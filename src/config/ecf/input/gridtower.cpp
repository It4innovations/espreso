
#include "gridtower.h"

#include "config/configuration.hpp"

espreso::GridTowerGeneratorConfiguration::GridTowerGeneratorConfiguration()
{
	direction = DIRECTION::X;
	composition = COMPOSITION::GLUED;

	REGISTER(direction, ECFMetaData()
			.setdescription({ "Direction of generated tower." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("X").setdescription("Grids are placed in x-direction."))
			.addoption(ECFOption().setname("Y").setdescription("Grids are placed in y-direction."))
			.addoption(ECFOption().setname("Z").setdescription("Grids are placed in z-direction.")));

	REGISTER(composition, ECFMetaData()
			.setdescription({ "Composition of generated cubes." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("GLUED").setdescription("Grids are glued into tower."))
			.addoption(ECFOption().setname("FREE").setdescription("Grids are kept untouched.")));

	REGISTER(grids, ECFMetaData()
			.setdescription({ "An index of grid in tower. Indices has to be continuous starting from 0.", "Description of grid." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
			.setpattern({ "0" }));
}


