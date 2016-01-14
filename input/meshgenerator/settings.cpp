
#include "settings.h"

using namespace esinput;

static std::vector<Description> createSetting()
{
	std::vector<Description> description(Settings::description);

	description.push_back({
		BOOLEAN_PARAMETER, "USE_METIS", "Use METIS for mesh partition."
	});

	return description;
};

std::vector<Description> Settings::description = createSetting();

Settings::Settings(int argc, char** argv,size_t index, size_t size)
: index(index), size(size)
{
	Configuration configuration(Settings::description, argc, argv);

	useMetis = configuration.value<eslocal>("USE_METIS", false);
}

Settings::Settings(size_t index, size_t size)
: index(index), size(size)
{
	useMetis = false;
}



