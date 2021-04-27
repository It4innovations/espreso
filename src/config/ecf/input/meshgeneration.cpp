
#include "meshgeneration.h"
#include "config/configuration.hpp"

using namespace espreso;

GMSHConfiguration::CharacteristicLength::CharacteristicLength()
{
	extern_from_boundary = 2;
	REGISTER(extern_from_boundary, ECFMetaData()
			.setdescription({ "Characteristic length external from boundary." })
			.setdatatype({ ECFDataType::FLOAT }));
	from_points = 0.9;
	REGISTER(from_points, ECFMetaData()
			.setdescription({ "Characteristic length from points." })
			.setdatatype({ ECFDataType::FLOAT }));
	from_curvature = 20;
	REGISTER(from_curvature, ECFMetaData()
			.setdescription({ "Characteristic length from curvature." })
			.setdatatype({ ECFDataType::FLOAT }));
	min = 2;
	REGISTER(min, ECFMetaData()
			.setdescription({ "Minimal characteristic length." })
			.setdatatype({ ECFDataType::FLOAT }));
	max = 2;
	REGISTER(max, ECFMetaData()
			.setdescription({ "Maximal characteristic length." })
			.setdatatype({ ECFDataType::FLOAT }));
}

GMSHConfiguration::GMSHConfiguration()
{
	algorithm3D = 7;
	REGISTER(algorithm3D, ECFMetaData()
			.setdescription({ "3D algorithm." })
			.setdatatype({ ECFDataType::INTEGER }));
	subdivisionAlgorithm = 1;
	REGISTER(subdivisionAlgorithm, ECFMetaData()
			.setdescription({ "Subdivision algorithm." })
			.setdatatype({ ECFDataType::INTEGER }));
	optimize = 0;
	REGISTER(optimize, ECFMetaData()
			.setdescription({ "optimize." })
			.setdatatype({ ECFDataType::INTEGER }));

	stl_angle = 40;
	REGISTER(stl_angle, ECFMetaData()
			.setdescription({ "Angle for STL meshing." })
			.setdatatype({ ECFDataType::INTEGER }));

	stl_precision = 1;
	REGISTER(stl_precision, ECFMetaData()
			.setdescription({ "Precision of STL meshing." })
			.setdatatype({ ECFDataType::INTEGER }));

	REGISTER(characteristic_length, ECFMetaData()
			.setdescription({ "Characteristic length." }));
}

MeshGenerationConfiguration::MeshGenerationConfiguration()
{
	REGISTER(gmsh_options, ECFMetaData()
			.setdescription({ "GMSH options." }));
}

