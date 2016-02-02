#include "settings.h"

using namespace esinput;

static std::vector<Description> createUniformSetting()
{
	std::vector<Description> description;
	std::vector<std::pair<std::string, std::string> > axis = {
			{"X", "x"},
			{"Y", "y"},
			{"Z", "z"}
	};

	for (size_t i = 0; i < axis.size(); i++) {
		description.push_back({
			INTEGER_PARAMETER, "SUBDOMAINS_" + axis[i].first, "Number of sub-domains in clusters in " + axis[i].second + "-axis."
		});
		description.push_back({
			INTEGER_PARAMETER, "ELEMENTS_" + axis[i].first, "Number of elements in clusters in " + axis[i].second + "-axis."
		});
		description.push_back({
			INTEGER_PARAMETER, "MATERIALS_" + axis[i].first, "Number of materials layers in " + axis[i].second + "-axis."
		});
	}

	description.push_back({ INTEGER_PARAMETER, "CORNER_COUNT", "The number of corners."});
	description.push_back({ BOOLEAN_PARAMETER, "CORNERS_IN_CORNERS", "Set corners in corner points."});
	description.push_back({ BOOLEAN_PARAMETER, "CORNERS_IN_EDGES", "Set corners on edges."});
	description.push_back({ BOOLEAN_PARAMETER, "CORNERS_IN_FACES", "Set corners on faces."});

	return description;
};

std::vector<Description> UniformSettings::description = createUniformSetting();

UniformSettings::UniformSettings(int argc, char** argv, size_t index, size_t size)
: Settings(index, size)
{
	Configuration configuration(UniformSettings::description, argc, argv);

	std::vector<std::string> axis = { "X", "Y", "Z" };
	for (size_t i = 0; i < axis.size(); i++) {
		subdomainsInCluster[i] = configuration.value<eslocal>("SUBDOMAINS_" + axis[i], 2);
		elementsInSubdomain[i] = configuration.value<eslocal>("ELEMENTS_" + axis[i], 5);
		materialsLayers[i] = configuration.value<eslocal>("MATERIALS_" + axis[i], 2);
	}

	cornerCount = configuration.value<eslocal>("CORNER_COUNT", 0);
	corners = configuration.value<bool>("CORNERS_IN_CORNERS", false);
	edges = configuration.value<bool>("CORNERS_IN_EDGES", false);
	faces = configuration.value<bool>("CORNERS_IN_FACES", false);
}

UniformSettings::UniformSettings(size_t index, size_t size)
: Settings(index, size)
{
	for (size_t i = 0; i < 3; i++) { // x, y, z
		subdomainsInCluster[i] = 2;
		elementsInSubdomain[i] = 5;
	}

	cornerCount = 0;
	corners = false;
	edges = false;
	faces = false;
}
