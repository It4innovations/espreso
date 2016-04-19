
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(CubeSettings &settings)
{
	for (size_t i = 0; i < 3; i++) {
		settings.clusters[i] = 1;
		settings.problemLength[i] = 30;
	}

	settings.boundaryCondition = std::vector<double>(6 * 2 * 3, std::numeric_limits<double>::infinity());
}

CubeSettings::CubeSettings(const Options &options, size_t index, size_t size)
: UniformSettings(options, index, size)
{
	defaultSettings(*this);
	ESINFO(OVERVIEW) << "Load cube setting from file " << options.path;

	std::vector<std::pair<std::string, std::string> > axis = {
			{ "X", "x" },
			{ "Y", "y" },
			{ "Z", "z" }
	};
	std::vector<std::pair<std::string, std::string> > properties = {
			{ "DIRICHLET", "Dirichlet" },
			{ "FORCES", "Force" }
	};
	std::vector<std::pair<std::string, std::string> > cube_faces = {
			{ "FRONT", "front" },
			{ "REAR", "rear" },
			{ "LEFT", "left" },
			{ "RIGHT", "right" },
			{ "TOP", "top" },
			{ "BOTTOM", "bottom" }
	};

	for (size_t i = 0; i < axis.size(); i++) {
		description.push_back({
			"CLUSTERS_" + axis[i].first, clusters[i], "Number of clusters in " + axis[i].second + "-axis."
		});
		description.push_back({
			"LENGTH_" + axis[i].first, problemLength[i], "Length of the cube in " + axis[i].second + "-axis."
		});
		for (size_t j = 0; j < properties.size(); j++) {
			for (size_t k = 0; k < cube_faces.size(); k++) {
				description.push_back({
					properties[j].first + "_" + cube_faces[k].first + "_" + axis[i].first, boundaryCondition[k * cube_faces.size() + j * properties.size() + i],
					properties[j].second + " on the " + cube_faces[k].second + " face in " + axis[i].second + "-axis."
				});
			}
		}
	}

	Configuration configuration(CubeSettings::description, options);
}

CubeSettings::CubeSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	defaultSettings(*this);
}


