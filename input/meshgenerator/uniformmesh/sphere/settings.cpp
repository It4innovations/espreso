
#include "settings.h"

using namespace espreso::input;

static std::vector<Description> createSphereSetting()
{
	std::vector<Description> description(UniformSettings::description);

	description.push_back({
		INTEGER_PARAMETER, "LAYERS", "Number of layers of the sphere."
	});
	description.push_back({
		DOUBLE_PARAMETER, "INNER_RADIUS", "Inner radius of the sphere."
	});
	description.push_back({
		DOUBLE_PARAMETER, "OUTER_RADIUS", "Outer radius of the sphere."
	});

	return description;
};

std::vector<Description> SphereSettings::description = createSphereSetting();

SphereSettings::SphereSettings(const Options &options, size_t index, size_t size)
: UniformSettings(options, index, size)
{
	ESINFO(OVERVIEW) << "Load sphere setting from file " << options.path;
	Configuration configuration(SphereSettings::description, options);

	layers = configuration.value<eslocal>("LAYERS", 1);
	innerRadius = configuration.value<double>("INNER_RADIUS", 9);
	outerRadius = configuration.value<double>("OUTER_RADIUS", 12);
}

SphereSettings::SphereSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	layers = 1;
	innerRadius = 9;
	outerRadius = 12;
}



