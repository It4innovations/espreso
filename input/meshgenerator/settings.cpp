
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(Settings &settings)
{
	settings.useMetis    = false;
	settings.shape       = CUBE;
	settings.elementType = HEXA8;
	settings.assembler   = LinearElasticity;

	settings.materials.push_back({7850, 2.1e11, 0.3});
	settings.materials.push_back({7850, 2.1e11, 0.3});
}

Settings::Settings(const Options &options, size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings(*this);
	description = {
		{ "USE_METIS"   , useMetis   , "Use METIS for mesh partition." },
		{ "SHAPE"       , shape      , "Generated shape. Supported values: 0 - CUBE, 1 - SPHERE" },
		{ "ELEMENT_TYPE", elementType, "The type of generated element. Supported values: <0, 7>" },

		{ "MAT1_DENSITY", materials[0].density     , "Density of the first material." },
		{ "MAT1_YOUNG"  , materials[0].youngModulus, "Young's modulus of the first material." },
		{ "MAT1_POISSON", materials[0].poissonRatio, "Poisson's ratio of the first material." },
		{ "MAT2_DENSITY", materials[1].density     , "Density of the first material." },
		{ "MAT2_YOUNG"  , materials[1].youngModulus, "Young's modulus of the first material." },
		{ "MAT2_POISSON", materials[1].poissonRatio, "Poisson's ratio of the first material." },

		{ "ASSEMBLER"   , assembler  , "Assembler type: 0 - LinearElasticity, 1 - Temperature" },
		{ "TIME_STEPS", config::assembler::timeSteps, "Number of time steps for transient problems."}
	};

	Configuration configuration(Settings::description, options);
}

Settings::Settings(size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings(*this);
}



