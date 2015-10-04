
#include "settings.h"

using namespace esinput;

std::vector<Description> SphereSettings::description = {
	{ INTEGER_PARAMETER, "SUBDOMAINS_X", "Number of sub-domains in clusters in x-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Y", "Number of sub-domains in clusters in y-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Z", "Number of sub-domains in clusters in z-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_X", "Number of elements in sub-domains in x-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Y", "Number of elements in sub-domains in y-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Z", "Number of elements in sub-domains in z-axis." }
};

SphereSettings::SphereSettings(int argc, char** argv)
{
	Configuration<SphereSettings> configuration(argc, argv);

	clusters[0] = clusters[1] = clusters[2] = 1;
	subdomainsInCluster[0] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	subdomainsInCluster[1] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	subdomainsInCluster[2] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	elementsInSubdomain[0] = configuration.value<eslocal>("ELEMENTS_X", 5);
	elementsInSubdomain[1] = configuration.value<eslocal>("ELEMENTS_X", 5);
	elementsInSubdomain[2] = configuration.value<eslocal>("ELEMENTS_X", 5);
	layers = configuration.value<eslocal>("LAYERS", 1);
}



