
#include "generator.h"

#include "config/configuration.hpp"

espreso::InputGeneratorConfiguration::InputGeneratorConfiguration()
{
	shape = INPUT_GENERATOR_SHAPE::GRID;

	REGISTER(shape, ECFMetaData()
			.setdescription({ "A generated shape." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("GRID").setdescription("A rectangular grid."))
			.addoption(ECFOption().setname("GRID_TOWER").setdescription("A rectangular grids in tower."))
			.addoption(ECFOption().setname("SPHERE").setdescription("A sphere.")));

	REGISTER(uniform_clusters, ECFMetaData()
			.setdescription({ "Turn ParMETIS decomposition on/off." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(uniform_domains, ECFMetaData()
			.setdescription({ "Turn METIS decomposition on/off." })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(grid, ECFMetaData()
			.setdescription({ "Settings of grid generator." })
			.allowonly([&] () { return shape == INPUT_GENERATOR_SHAPE::GRID; }));
	REGISTER(grid_tower, ECFMetaData()
			.setdescription({ "Settings of grid tower generator." })
			.allowonly([&] () { return shape == INPUT_GENERATOR_SHAPE::GRID_TOWER; }));
	REGISTER(sphere, ECFMetaData()
			.setdescription({ "Settings of sphere generator." })
			.allowonly([&] () { return shape == INPUT_GENERATOR_SHAPE::SPHERE; }));

	uniform_clusters = uniform_domains = true;
}



