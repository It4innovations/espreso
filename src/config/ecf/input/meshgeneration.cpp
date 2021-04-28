
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

NGLibConfiguration::NGLibConfiguration()
{
	uselocalh = 1;
	REGISTER(uselocalh, ECFMetaData()
			.setdescription({ "Switch to enable / disable usage of local mesh size modifiers." })
			.setdatatype({ ECFDataType::INTEGER }));
	maxh = 1000;
	REGISTER(maxh, ECFMetaData()
			.setdescription({ "Maximum global mesh size allowed." })
			.setdatatype({ ECFDataType::FLOAT }));
	minh = 0;
	REGISTER(minh, ECFMetaData()
			.setdescription({ "Minimum global mesh size allowed." })
			.setdatatype({ ECFDataType::FLOAT }));

	fineness = .5;
	REGISTER(fineness, ECFMetaData()
			.setdescription({ "Mesh density: 0...1 (0 => coarse; 1 => fine)." })
			.setdatatype({ ECFDataType::FLOAT }));
	grading = .3;
	REGISTER(grading, ECFMetaData()
			.setdescription({ "Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)." })
			.setdatatype({ ECFDataType::FLOAT }));
	elementsperedge = 2;
	REGISTER(elementsperedge, ECFMetaData()
			.setdescription({ "Number of elements to generate per edge of the geometry." })
			.setdatatype({ ECFDataType::FLOAT }));
	elementspercurve = 2;
	REGISTER(elementspercurve, ECFMetaData()
			.setdescription({ "Elements to generate per curvature radius." })
			.setdatatype({ ECFDataType::FLOAT }));

	closeedgeenable = 0;
	REGISTER(closeedgeenable, ECFMetaData()
			.setdescription({ "Enable / Disable mesh refinement at close edges." })
			.setdatatype({ ECFDataType::INTEGER }));
	closeedgefact = 2;
	REGISTER(closeedgefact, ECFMetaData()
			.setdescription({ "Factor to use for refinement at close edges (larger => finer)." })
			.setdatatype({ ECFDataType::FLOAT }));

	minedgelenenable = 0;
	REGISTER(minedgelenenable, ECFMetaData()
			.setdescription({ "Enable / Disable user defined minimum edge length for edge subdivision." })
			.setdatatype({ ECFDataType::INTEGER }));
	minedgelen = 0;
	REGISTER(minedgelen, ECFMetaData()
			.setdescription({ "Minimum edge length to use while subdividing the edges (default = 1e-4)." })
			.setdatatype({ ECFDataType::FLOAT }));

	second_order = 0;
	REGISTER(second_order, ECFMetaData()
			.setdescription({ "Generate second-order surface and volume elements." })
			.setdatatype({ ECFDataType::INTEGER }));
	quad_dominated = 0;
	REGISTER(quad_dominated, ECFMetaData()
			.setdescription({ "Creates a Quad-dominated mesh." })
			.setdatatype({ ECFDataType::INTEGER }));
	optsurfmeshenable = 1;
	REGISTER(optsurfmeshenable, ECFMetaData()
			.setdescription({ "Enable / Disable automatic surface mesh optimization." })
			.setdatatype({ ECFDataType::INTEGER }));
	optvolmeshenable = 1;
	REGISTER(optvolmeshenable, ECFMetaData()
			.setdescription({ "Enable / Disable automatic volume mesh optimization." })
			.setdatatype({ ECFDataType::INTEGER }));
	optsteps_3d = 3;
	REGISTER(optsteps_3d, ECFMetaData()
			.setdescription({ "Number of optimize steps to use for 3-D mesh optimization." })
			.setdatatype({ ECFDataType::INTEGER }));
	optsteps_2d = 2;
	REGISTER(optsteps_2d, ECFMetaData()
			.setdescription({ "Number of optimize steps to use for 2-D mesh optimization." })
			.setdatatype({ ECFDataType::INTEGER }));
	invert_tets = 0;
	REGISTER(invert_tets, ECFMetaData()
			.setdescription({ "Invert all the volume elements." })
			.setdatatype({ ECFDataType::INTEGER }));
	invert_trigs = 0;
	REGISTER(invert_trigs, ECFMetaData()
			.setdescription({ "Invert all the surface triangle elements." })
			.setdatatype({ ECFDataType::INTEGER }));
	check_overlap = 1;
	REGISTER(check_overlap, ECFMetaData()
			.setdescription({ "Check for overlapping surfaces during Surface meshing." })
			.setdatatype({ ECFDataType::INTEGER }));
	check_overlapping_boundary = 1;
	REGISTER(check_overlapping_boundary, ECFMetaData()
			.setdescription({ "Check for overlapping surface elements before volume meshing." })
			.setdatatype({ ECFDataType::INTEGER }));
}

MeshGenerationConfiguration::MeshGenerationConfiguration()
{
	REGISTER(gmsh_options, ECFMetaData()
			.setdescription({ "GMSH options." }));

	REGISTER(nglib_options, ECFMetaData()
			.setdescription({ "NGLib options." }));
}

