
#include "noderegion.h"
#include "config/configuration.hpp"

using namespace espreso;

InputNodeRegionConfiguration::InputNodeRegionConfiguration()
{
	shape = SHAPE::CIRCLE;
	REGISTER(shape, ECFMetaData()
			.setdescription({ "Region shape." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CIRCLE").setdescription("Region from a circle."))
			.addoption(ECFOption().setname("BLOCK").setdescription("Region from a block.")));

	x = 0;
	REGISTER(x, ECFMetaData()
			.setdescription({ "X" })
			.setdatatype({ ECFDataType::FLOAT }));
	y = 0;
	REGISTER(y, ECFMetaData()
			.setdescription({ "Y" })
			.setdatatype({ ECFDataType::FLOAT }));
	z = 0;
	REGISTER(z, ECFMetaData()
			.setdescription({ "Z" })
			.setdatatype({ ECFDataType::FLOAT }));

	lx = 0;
	REGISTER(lx, ECFMetaData()
			.setdescription({ "LENGTH X" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return shape == SHAPE::BLOCK; }));
	ly = 0;
	REGISTER(ly, ECFMetaData()
			.setdescription({ "LENGTH Y" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return shape == SHAPE::BLOCK; }));
	lz = 0;
	REGISTER(lz, ECFMetaData()
			.setdescription({ "LENGTH Z" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return shape == SHAPE::BLOCK; }));

	radius = 0;
	REGISTER(radius, ECFMetaData()
			.setdescription({ "RADIUS" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return shape == SHAPE::CIRCLE; }));
}



