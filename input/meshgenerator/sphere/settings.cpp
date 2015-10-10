
#include "settings.h"

using namespace esinput;

std::vector<Description> SphereSettings::description = {
	{ INTEGER_PARAMETER, "SUBDOMAINS_X", "Number of sub-domains in clusters in x-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Y", "Number of sub-domains in clusters in y-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Z", "Number of sub-domains in clusters in z-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_X", "Number of elements in sub-domains in x-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Y", "Number of elements in sub-domains in y-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Z", "Number of elements in sub-domains in z-axis." },
	{ INTEGER_PARAMETER, "LAYERS", "Number of layers of the sphere." },
	{ DOUBLE_PARAMETER, "INNER_RADIUS", "Inner radius of the sphere." },
	{ DOUBLE_PARAMETER, "OUTER_RADIUS", "Outer radius of the sphere." },

	{ INTEGER_PARAMETER, "CORNER_COUNT", "The number of corners."},
	{ BOOLEAN_PARAMETER, "CORNERS_IN_CORNERS", "Set corners in corner points."},
	{ BOOLEAN_PARAMETER, "CORNERS_IN_EDGES", "Set corners on edges."},
	{ BOOLEAN_PARAMETER, "CORNERS_IN_FACES", "Set corners on faces."}

};

SphereSettings::SphereSettings(int argc, char** argv)
{
	Configuration<SphereSettings> configuration(argc, argv);

	subdomainsInCluster[0] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	subdomainsInCluster[1] = configuration.value<eslocal>("SUBDOMAINS_Y", 2);
	subdomainsInCluster[2] = configuration.value<eslocal>("SUBDOMAINS_Z", 2);
	elementsInSubdomain[0] = configuration.value<eslocal>("ELEMENTS_X", 5);
	elementsInSubdomain[1] = configuration.value<eslocal>("ELEMENTS_Y", 5);
	elementsInSubdomain[2] = configuration.value<eslocal>("ELEMENTS_Z", 5);
	layers = configuration.value<eslocal>("LAYERS", 1);
	innerRadius = configuration.value<double>("INNER_RADIUS", 9);
	outerRadius = configuration.value<double>("OUTER_RADIUS", 12);

	cornerCount = configuration.value<eslocal>("CORNER_COUNT", 0);
	corners = configuration.value<bool>("CORNERS_IN_CORNERS", true);
	edges = configuration.value<bool>("CORNERS_IN_EDGES", false);
	faces = configuration.value<bool>("CORNERS_IN_FACES", false);
}



