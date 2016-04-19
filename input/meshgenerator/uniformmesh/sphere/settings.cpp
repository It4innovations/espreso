
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(SphereSettings &settings)
{
	settings.layers = 1;
	settings.grid = 1;
	settings.innerRadius = 9;
	settings.outerRadius = 12;

	settings.boundaryCondition = std::vector<double>(2 * 2 * 3);
}

SphereSettings::SphereSettings(const Options &options, size_t index, size_t size)
: UniformSettings(options, index, size)
{
	ESINFO(OVERVIEW) << "Load sphere setting from file " << options.path;

	defaultSettings(*this);

	description.push_back({
		"LAYERS", layers, "Number of layers of the sphere."
	});
	description.push_back({
		"GRID", grid, "Grid size of one side of the sphere."
	});
	description.push_back({
		"INNER_RADIUS", innerRadius, "Inner radius of the sphere."
	});
	description.push_back({
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
				description.push_back({
					properties[j].first + "_" + sphere_faces[k].first + "_" + axis[i].first, boundaryCondition[k * sphere_faces.size() + j * properties.size() + i],
					properties[j].second + " on the " + sphere_faces[k].second + " face in " + axis[i].second + "-axis."
				});
			}
		}
	}

	Configuration configuration(SphereSettings::description, options);
}

SphereSettings::SphereSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	defaultSettings(*this);
}



