
#ifndef SRC_CONFIGURATION_INPUT_INPUTGENERATORGRID_H_
#define SRC_CONFIGURATION_INPUT_INPUTGENERATORGRID_H_

#include "inputgeneratorelements.h"

namespace espreso {

struct GridConfiguration: public Configuration {

	OPTION(ELEMENT_TYPE, element_type, "Type of generated element", ELEMENT_TYPE::HEXA8, OPTIONS({
		{ "HEXA8"    , ELEMENT_TYPE::HEXA8    , "Hexahedron."},
		{ "HEXA20"   , ELEMENT_TYPE::HEXA20   , "Hexahedron with midpoints."},
		{ "TETRA4"   , ELEMENT_TYPE::TETRA4   , "Tetrahedron."},
		{ "TETRA10"  , ELEMENT_TYPE::TETRA10  , "Tetrahedron with midpoints."},
		{ "PRISMA6"  , ELEMENT_TYPE::PRISMA6  , "Prisma."},
		{ "PRISMA15" , ELEMENT_TYPE::PRISMA15 , "Prisma with midpoints."},
		{ "PYRAMID5" , ELEMENT_TYPE::PYRAMID5 , "Pyramid."},
		{ "PYRAMID13", ELEMENT_TYPE::PYRAMID13, "Pyramid with midpoints."},

		{ "SQUARE4"  , ELEMENT_TYPE::SQUARE4  , "Square."},
		{ "SQUARE8"  , ELEMENT_TYPE::SQUARE8  , "Square with midpoints."},
		{ "TRIANGLE3", ELEMENT_TYPE::TRIANGLE3, "Triangle."},
		{ "TRIANGLE6", ELEMENT_TYPE::TRIANGLE6, "Triangle with midpoints."},
	}));

	PARAMETER(double, start_x, "x-coordinate of grid starting point.", 0);
	PARAMETER(double, start_y, "y-coordinate of grid starting point.", 0);
	PARAMETER(double, start_z, "z-coordinate of grid starting point.", 0);
	PARAMETER(double, length_x, "x-length of generated grid.", 1);
	PARAMETER(double, length_y, "y-length of generated grid.", 1);
	PARAMETER(double, length_z, "z-length of generated grid.", 1);

	PARAMETER(std::string, projection_x, "Projection of x-coordinate.", "x");
	PARAMETER(std::string, projection_y, "Projection of y-coordinate.", "y");
	PARAMETER(std::string, projection_z, "Projection of z-coordinate.", "z");
	PARAMETER(std::string, rotation_x, "Rotation of x-coordinate.", "0");
	PARAMETER(std::string, rotation_y, "Rotation of y-coordinate.", "0");
	PARAMETER(std::string, rotation_z, "Rotation of z-coordinate.", "0");


	PARAMETER(double, blocks_x, "Number of blocks in x-direction of a grid.", 1);
	PARAMETER(double, blocks_y, "Number of blocks in y-direction of a grid.", 1);
	PARAMETER(double, blocks_z, "Number of blocks in z-direction of a grid.", 1);

	PARAMETER(double, clusters_x, "Number of clusters in x-direction of each grid square.", 1);
	PARAMETER(double, clusters_y, "Number of clusters in y-direction of each grid square.", 1);
	PARAMETER(double, clusters_z, "Number of clusters in z-direction of each grid square.", 1);

	PARAMETER(double, domains_x, "Number of domains in x-direction of each cluster.", 2);
	PARAMETER(double, domains_y, "Number of domains in y-direction of each cluster.", 2);
	PARAMETER(double, domains_z, "Number of domains in z-direction of each cluster.", 2);

	PARAMETER(double, elements_x, "Number of elements in x-direction of each domain.", 5);
	PARAMETER(double, elements_y, "Number of elements in y-direction of each domain.", 5);
	PARAMETER(double, elements_z, "Number of elements in z-direction of each domain.", 5);

	PARAMETER(bool, uniform_decomposition, "Grid is uniformly decomposed", true);

	SUBMAP(size_t, bool, blocks, "List of grid blocks [<INDEX> <VALUE>]. Where value indicate if a block will be generated.", "0", true);

	SUBMAP(std::string, std::string, nodes, "List of nodes regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, edges, "List of edges regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, faces, "List of faces regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, elements, "List of elements regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL;CHESSBOARD_WHITE;CHESSBOARD_BLACK;NOT_SELECTED]>");

	PARAMETER(size_t, chessboard_size, "Number of squares of chessboard in one direction", 2);
};

}



#endif /* SRC_CONFIGURATION_INPUT_INPUTGENERATORGRID_H_ */
