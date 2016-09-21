
#include "settings.h"

using namespace espreso::input;

void SphereSettings::defaultSphereSettings()
{
	layers = 1;
	grid = 1;
	innerRadius = 9;
	outerRadius = 12;

	parameters.push_back({
		prefix + "LAYERS", layers, "Number of layers of the sphere."
	});
	parameters.push_back({
		prefix + "GRID", grid, "Grid size of one side of the sphere."
	});
	parameters.push_back({
		prefix + "INNER_RADIUS", innerRadius, "Inner radius of the sphere."
	});
	parameters.push_back({
		prefix + "OUTER_RADIUS", outerRadius, "Outer radius of the sphere."
	});
}

SphereSettings::SphereSettings(const Configuration &configuration, size_t index, size_t size, std::string prefix)
: UniformSettings(index, size, prefix)
{
	ESINFO(OVERVIEW) << "Load sphere setting from file " << configuration.path;

	defaultSphereSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

SphereSettings::SphereSettings(size_t index, size_t size, std::string prefix)
: UniformSettings(index, size, prefix)
{
	defaultSphereSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
}



