
#include "settings.h"

using namespace espreso::input;

static std::vector<Description> createSetting()
{
	std::vector<Description> description(Settings::description);

	description.push_back({
		BOOLEAN_PARAMETER, "USE_METIS", "Use METIS for mesh partition."
	});

	description.push_back({
		INTEGER_PARAMETER, "SHAPE", "Generated shape. Supported values: 0 - CUBE, 1 - SPHERE"
	});
	description.push_back({
		INTEGER_PARAMETER, "ELEMENT_TYPE", "The type of generated element. Supported values: <0, 7>"
	});

	description.push_back({
		INTEGER_PARAMETER, "ASSEMBLER", "Assembler type: 0 - LinearElasticity, 1 - Temperature"
	});

	return description;
};

std::vector<Description> Settings::description = createSetting();

Settings::Settings(const Options &options, size_t index, size_t size)
: index(index), size(size)
{
	Configuration configuration(Settings::description, options);

	useMetis = configuration.value("USE_METIS", false);
	shape = configuration.value("SHAPE", 0);
	elementType = configuration.value("ELEMENT_TYPE", 0);
	assembler = configuration.value("ASSEMBLER", 0);
}

Settings::Settings(size_t index, size_t size)
: index(index), size(size)
{
	useMetis = false;
	shape = 0;
	elementType = 0;
	assembler = 0;
}



