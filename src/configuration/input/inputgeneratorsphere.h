
#ifndef SRC_CONFIGURATION_INPUT_INPUTGENERATORSPHERE_H_
#define SRC_CONFIGURATION_INPUT_INPUTGENERATORSPHERE_H_

#include "inputgeneratorelements.h"

namespace espreso {

struct SphereConfiguration: public Configuration {

	OPTION(ELEMENT_TYPE, element_type, "Type of generated element", ELEMENT_TYPE::HEXA8, OPTIONS({
		{ "HEXA8"    , ELEMENT_TYPE::HEXA8    , "Hexahedron."},
		{ "HEXA20"   , ELEMENT_TYPE::HEXA20   , "Hexahedron with midpoints."},
		{ "TETRA4"   , ELEMENT_TYPE::TETRA4   , "Tetrahedron."},
		{ "TETRA10"  , ELEMENT_TYPE::TETRA10  , "Tetrahedron with midpoints."},
		{ "PRISMA6"  , ELEMENT_TYPE::PRISMA6  , "Prisma."},
		{ "PRISMA15" , ELEMENT_TYPE::PRISMA15 , "Prisma with midpoints."},
		{ "PYRAMID5" , ELEMENT_TYPE::PYRAMID5 , "Pyramid."},
		{ "PYRAMID13", ELEMENT_TYPE::PYRAMID13, "Pyramid with midpoints."},
	}));

	PARAMETER(double, inner_radius, "Inner radius of generated sphere.", 5);
	PARAMETER(double, outer_radius, "Outer radius of generated sphere.", 10);

	PARAMETER(double, clusters, "Number of clusters in x,y-directions of each square.", 1);
	PARAMETER(double, layers, "Number of clusters in z-direction.", 1);

	PARAMETER(double, domains_x, "Number of domains in x-direction of each cluster.", 2);
	PARAMETER(double, domains_y, "Number of domains in y-direction of each cluster.", 2);
	PARAMETER(double, domains_z, "Number of domains in z-direction of each cluster.", 2);

	PARAMETER(double, elements_x, "Number of elements in x-direction of each domain.", 5);
	PARAMETER(double, elements_y, "Number of elements in y-direction of each domain.", 5);
	PARAMETER(double, elements_z, "Number of elements in z-direction of each domain.", 5);

	PARAMETER(bool, uniform_decomposition, "Grid is uniformly decomposed", true);

	SUBMAP(std::string, std::string, nodes, "List of nodes regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, edges, "List of edges regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, faces, "List of faces regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, elements, "List of elements regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
};

}



#endif /* SRC_CONFIGURATION_INPUT_INPUTGENERATORSPHERE_H_ */
