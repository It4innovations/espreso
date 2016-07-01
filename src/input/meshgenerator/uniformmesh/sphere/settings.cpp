
#include "settings.h"

using namespace espreso::input;

void SphereSettings::defaultSphereSettings()
{
	layers = 1;
	grid = 1;
	innerRadius = 9;
	outerRadius = 12;

	boundaryCondition = std::vector<double>(2 * 2 * 3, std::numeric_limits<double>::infinity());

	parameters.push_back({
		"LAYERS", layers, "Number of layers of the sphere."
	});
	parameters.push_back({
		"GRID", grid, "Grid size of one side of the sphere."
	});
	parameters.push_back({
		"INNER_RADIUS", innerRadius, "Inner radius of the sphere."
	});
	parameters.push_back({
		"OUTER_RADIUS", outerRadius, "Outer radius of the sphere."
	});

	std::vector<std::pair<std::string, std::string> > axis = {
			{ "X", "x" },
			{ "Y", "y" },
			{ "Z", "z" }
	};
	std::vector<std::pair<std::string, std::string> > properties = {
			{ "DIRICHLET", "Dirichlet" },
			{ "FORCES", "Force" }
	};
	std::vector<std::pair<std::string, std::string> > sphere_faces = {
			{ "INNER", "inner" },
			{ "OUTER", "outer" }
	};

	for (size_t i = 0; i < axis.size(); i++) {
		for (size_t j = 0; j < properties.size(); j++) {
			for (size_t k = 0; k < sphere_faces.size(); k++) {
				parameters.push_back({
					properties[j].first + "_" + sphere_faces[k].first + "_" + axis[i].first, boundaryCondition[k * properties.size() * axis.size() + j * properties.size() + i],
					properties[j].second + " on the " + sphere_faces[k].second + " face in " + axis[i].second + "-axis."
				});
			}
		}
	}
}

SphereSettings::SphereSettings(const Configuration &configuration, size_t index, size_t size)
: UniformSettings(index, size)
{
	ESINFO(OVERVIEW) << "Load sphere setting from file " << configuration.path;

	defaultSphereSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

SphereSettings::SphereSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	defaultSphereSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
}



