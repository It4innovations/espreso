
#include "settings.h"

using namespace esinput;

std::vector<Description> FactorySettings::description = {
	{ INTEGER_PARAMETER, "SHAPE", "Generated shape. Supported values: 0 - CUBE, 1 - SPHERE" },
	{ INTEGER_PARAMETER, "ELEMENT_TYPE", "The type of generated element. Supported values: <0, 7>" }
};

FactorySettings::FactorySettings(int argc, char** argv)
{
	Configuration<CubeSettings> configuration(argc, argv);

	shape = configuration.value<eslocal>("SHAPE", 0);
	elementType = configuration.value<eslocal>("ELEMENT_TYPE", 0);
}

