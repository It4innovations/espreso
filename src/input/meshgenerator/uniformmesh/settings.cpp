#include "settings.h"

using namespace espreso::input;

void UniformSettings::defaultUniformSettings()
{
	for (size_t i = 0; i < 3; i++) { // x, y, z
		subdomainsInCluster[i] = 2;
		elementsInSubdomain[i] = 5;
		materialsLayers[i] = 1;
	}

	cornerCount = 0;
	corners     = true;
	edges       = false;
	faces       = false;

	parameters = {
		{ prefix + "SUBDOMAINS_X", subdomainsInCluster[0], "Number of sub-domains in a cluster in x-axis."},
		{ prefix + "SUBDOMAINS_Y", subdomainsInCluster[1], "Number of sub-domains in a cluster in y-axis."},
		{ prefix + "SUBDOMAINS_Z", subdomainsInCluster[2], "Number of sub-domains in a cluster in z-axis."},

		{ prefix + "ELEMENTS_X", elementsInSubdomain[0], "Number of elements in a sub-domain in x-axis."},
		{ prefix + "ELEMENTS_Y", elementsInSubdomain[1], "Number of elements in a sub-domain in y-axis."},
		{ prefix + "ELEMENTS_Z", elementsInSubdomain[2], "Number of elements in a sub-domain in z-axis."},

		{ prefix + "MATERIALS_X", materialsLayers[0], "Number of materials layers in x-axis."},
		{ prefix + "MATERIALS_Y", materialsLayers[1], "Number of materials layers in y-axis."},
		{ prefix + "MATERIALS_Z", materialsLayers[2], "Number of materials layers in z-axis."},

		{ prefix + "CORNER_COUNT"      , cornerCount, "The number of corners."},
		{ prefix + "CORNERS_IN_CORNERS", corners    , "Set corners to corners points."},
		{ prefix + "CORNERS_IN_EDGES"  , edges      , "Set corners on edges."},
		{ prefix + "CORNERS_IN_FACES"  , faces      , "Set corners on faces."},
	};
}

UniformSettings::UniformSettings(const Configuration &configuration, size_t index, size_t size, std::string prefix)
: Settings(index, size, prefix)
{
	defaultUniformSettings();
	parameters.insert(parameters.end(), Settings::parameters.begin(), Settings::parameters.end());
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

UniformSettings::UniformSettings(size_t index, size_t size, std::string prefix)
: Settings(index, size, prefix)
{
	defaultUniformSettings();
	parameters.insert(parameters.end(), Settings::parameters.begin(), Settings::parameters.end());
}
