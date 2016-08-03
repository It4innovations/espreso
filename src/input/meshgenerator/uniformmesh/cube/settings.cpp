
#include "settings.h"

using namespace espreso::input;

void CubeSettings::defaultCubeSettings()
{
	for (size_t i = 0; i < 3; i++) {
		clusters[i] = 1;
		problemLength[i] = 30;
	}

	boundaryCondition = std::vector<double>(6 * 2 * 3, std::numeric_limits<double>::infinity());

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
		parameters.push_back({
			"CLUSTERS_" + axis[i].first, clusters[i], "Number of clusters in " + axis[i].second + "-axis."
		});
		parameters.push_back({
			"LENGTH_" + axis[i].first, problemLength[i], "Length of the cube in " + axis[i].second + "-axis."
		});
		for (size_t j = 0; j < properties.size(); j++) {
			for (size_t k = 0; k < cube_faces.size(); k++) {
				parameters.push_back({
					properties[j].first + "_" + cube_faces[k].first + "_" + axis[i].first, boundaryCondition[k * properties.size() * axis.size() + j * properties.size() + i],
					properties[j].second + " on the " + cube_faces[k].second + " face in " + axis[i].second + "-axis."
				});
			}
		}
	}
}

CubeSettings::CubeSettings(const Configuration &configuration, size_t index, size_t size)
: UniformSettings(index, size)
{
	defaultCubeSettings();
	ESINFO(OVERVIEW) << "Load cube setting from file " << configuration.path;
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

CubeSettings::CubeSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	defaultCubeSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
}


