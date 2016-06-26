
#include "settings.h"

using namespace espreso::input;

void PlaneSettings::defaultPlaneSettings()
{
	for (size_t i = 0; i < 2; i++) {
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
	std::vector<std::pair<std::string, std::string> > plane_faces = {
			{ "FRONT", "not supported" },
			{ "REAR", "not supported" },
			{ "LEFT", "left" },
			{ "RIGHT", "right" },
			{ "TOP", "top" },
			{ "BOTTOM", "bottom" }
	};

	for (size_t i = 0; i < 2; i++) {
		parameters.push_back({
			"CLUSTERS_" + axis[i].first, clusters[i], "Number of clusters in " + axis[i].second + "-axis."
		});
		parameters.push_back({
			"LENGTH_" + axis[i].first, problemLength[i], "Length of the cube in " + axis[i].second + "-axis."
		});
	}

	for (size_t i = 0; i < axis.size(); i++) {
		for (size_t j = 0; j < properties.size(); j++) {
			for (size_t k = 0; k < plane_faces.size(); k++) {
				parameters.push_back({
					properties[j].first + "_" + plane_faces[k].first + "_" + axis[i].first, boundaryCondition[k * properties.size() * axis.size() + j * properties.size() + i],
					properties[j].second + " on the " + plane_faces[k].second + " face in " + axis[i].second + "-axis."
				});
			}
		}
	}
}

PlaneSettings::PlaneSettings(const Configuration &configuration, size_t index, size_t size)
: CubeSettings(index, size)
{
	parameters.clear();
	defaultPlaneSettings();
	ESINFO(OVERVIEW) << "Load plane setting from file " << configuration.path;
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	ParametersReader::configuration(configuration, parameters);
	subdomainsInCluster[2] = 1;
	elementsInSubdomain[2] = 1;
}

PlaneSettings::PlaneSettings(size_t index, size_t size)
: CubeSettings(index, size)
{
	parameters.clear();
	defaultPlaneSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	subdomainsInCluster[2] = 1;
	elementsInSubdomain[2] = 1;
}


