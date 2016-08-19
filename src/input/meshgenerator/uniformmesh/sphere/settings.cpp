
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(SphereSettings &settings)
{
	settings.layers = 1;
	settings.grid = 1;
	settings.innerRadius = 9;
	settings.outerRadius = 12;

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

	parameters.push_back({ prefix + "DIRICHLET", dirichlet, "Dirichlet boundary conditions" });
	parameters.push_back({ prefix + "FORCES", forces, "Boundary forces" });
}

SphereSettings::SphereSettings(const Configuration &configuration, size_t index, size_t size, std::string prefix)
: UniformSettings(index, size, prefix)
{
	ESINFO(OVERVIEW) << "Load sphere setting from file " << configuration.path;

	Configuration configuration(SphereSettings::description, options);
}

SphereSettings::SphereSettings(size_t index, size_t size, std::string prefix)
: UniformSettings(index, size, prefix)
{
	defaultSettings(*this);
}



