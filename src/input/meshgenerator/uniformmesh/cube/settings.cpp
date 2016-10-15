
#include "settings.h"

using namespace espreso::input;

void CubeSettings::defaultCubeSettings()
{
	for (size_t i = 0; i < 3; i++) {
		clusters[i] = 1;
		problemOrigin[i] = 0;
		problemLength[i] = 30;
	}

	std::vector<std::pair<std::string, std::string> > axis = {
			{ "X", "x" },
			{ "Y", "y" },
			{ "Z", "z" }
	};

	size_t verbosity = 1;

	for (size_t i = 0; i < axis.size(); i++) {
		parameters.push_back({
			prefix + "CLUSTERS_" + axis[i].first, clusters[i], "Number of clusters in " + axis[i].second + "-axis.", verbosity
		});
	}
	for (size_t i = 0; i < axis.size(); i++) {
		parameters.push_back({
			prefix + "ORIGIN_" + axis[i].first, problemOrigin[i], "Length of the cube in " + axis[i].second + "-axis.", verbosity
		});
	}
	for (size_t i = 0; i < axis.size(); i++) {
		parameters.push_back({
			prefix + "LENGTH_" + axis[i].first, problemLength[i], "Length of the cube in " + axis[i].second + "-axis.", verbosity
		});
	}

	for (size_t i = 0; i < axis.size(); i++) {
		parameters.push_back({
			prefix + "PROJECTION_" + axis[i].first, projections[i], "Projection of " + axis[i].second + "-axis.", verbosity
		});
	}
	for (size_t i = 0; i < axis.size(); i++) {
		parameters.push_back({
			prefix + "ROTATION_" + axis[i].first, rotations[i], "Projection of " + axis[i].second + "-axis.", verbosity
		});
	}
}

CubeSettings::CubeSettings(const ArgsConfiguration &configuration, size_t index, size_t size, std::string prefix)
: UniformSettings(index, size, prefix)
{
	defaultCubeSettings();
	ESINFO(OVERVIEW) << "Load cube setting from file " << configuration.path;
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

CubeSettings::CubeSettings(size_t index, size_t size, std::string prefix)
: UniformSettings(index, size, prefix)
{
	defaultCubeSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
}


