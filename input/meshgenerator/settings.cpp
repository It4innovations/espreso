
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(Settings &settings)
{
	settings.useMetis    = false;
	settings.shape       = CUBE;
	settings.elementType = HEXA8;
	settings.assembler   = LinearElasticity;
}

Settings::Settings(const Options &options, size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings(*this);
	description = {
		{ "USE_METIS"   , useMetis   , "Use METIS for mesh partition." },
		{ "SHAPE"       , shape      , "Generated shape. Supported values: 0 - CUBE, 1 - SPHERE" },
		{ "ELEMENT_TYPE", elementType, "The type of generated element. Supported values: <0, 7>" },
		{ "ASSEMBLER"   , assembler  , "Assembler type: 0 - LinearElasticity, 1 - Temperature" }
	};

	Configuration configuration(Settings::description, options);
}

Settings::Settings(size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings(*this);
}



