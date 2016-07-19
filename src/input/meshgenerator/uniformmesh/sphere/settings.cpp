
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(SphereSettings &settings)
{
	settings.layers = 1;
	settings.grid = 1;
	settings.innerRadius = 9;
	settings.outerRadius = 12;

	parameters.push_back({
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

	parameters.push_back({ "DIRICHLET", dirichlet, "Dirichlet boundary conditions" });
	parameters.push_back({ "FORCES", forces, "Boundary forces" });
}

SphereSettings::SphereSettings(const Configuration &configuration, size_t index, size_t size)
: UniformSettings(index, size)
{
	ESINFO(OVERVIEW) << "Load sphere setting from file " << configuration.path;

	Configuration configuration(SphereSettings::description, options);
}

SphereSettings::SphereSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	defaultSettings(*this);
}



