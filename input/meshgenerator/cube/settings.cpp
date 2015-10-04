
#include "settings.h"

using namespace esinput;

std::vector<Description> CubeSettings::description = {
	{ INTEGER_PARAMETER, "CLUSTERS_X", "Number of clusters in x-axis." },
	{ INTEGER_PARAMETER, "CLUSTERS_Y", "Number of clusters in y-axis." },
	{ INTEGER_PARAMETER, "CLUSTERS_Z", "Number of clusters in z-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_X", "Number of sub-domains in clusters in x-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Y", "Number of sub-domains in clusters in y-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Z", "Number of sub-domains in clusters in z-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_X", "Number of elements in sub-domains in x-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Y", "Number of elements in sub-domains in y-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Z", "Number of elements in sub-domains in z-axis." },
	{ DOUBLE_PARAMETER, "LENGHT_X", "Length of the cube in x-axis." },
	{ DOUBLE_PARAMETER, "LENGHT_Y", "Length of the cube in y-axis." },
	{ DOUBLE_PARAMETER, "LENGHT_Z", "Length of the cube in z-axis." }
};

CubeSettings::CubeSettings(int argc, char** argv)
{
	Configuration<CubeSettings> configuration(argc, argv);

	clusters[0] = configuration.value<eslocal>("CLUSTERS_X", 1);
	clusters[1] = configuration.value<eslocal>("CLUSTERS_Y", 1);
	clusters[2] = configuration.value<eslocal>("CLUSTERS_Z", 1);
	subdomainsInCluster[0] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	subdomainsInCluster[1] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	subdomainsInCluster[2] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	elementsInSubdomain[0] = configuration.value<eslocal>("ELEMENTS_X", 5);
	elementsInSubdomain[1] = configuration.value<eslocal>("ELEMENTS_X", 5);
	elementsInSubdomain[2] = configuration.value<eslocal>("ELEMENTS_X", 5);
	problemLength[0] = configuration.value<double>("LENGTH_X", 20);
	problemLength[1] = configuration.value<double>("LENGTH_Y", 20);
	problemLength[2] = configuration.value<double>("LENGTH_Z", 20);
}


