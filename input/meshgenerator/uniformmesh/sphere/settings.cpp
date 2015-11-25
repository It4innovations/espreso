
#include "settings.h"

using namespace esinput;

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

SphereSettings::SphereSettings(int argc, char** argv,size_t index, size_t size)
: UniformSettings(argc, argv, index, size)
{
	Configuration configuration(SphereSettings::description, argc, argv);

	layers = configuration.value<eslocal>("LAYERS", 1);
	innerRadius = configuration.value<double>("INNER_RADIUS", 9);
	outerRadius = configuration.value<double>("OUTER_RADIUS", 12);
}



