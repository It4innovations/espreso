
#include "settings.h"

using namespace espreso::input;

void Settings::defaultSettings()
{
	useMetis    = false;
	shape       = CUBE;
	elementType = HEXA8;
	assembler   = LinearElasticity;

	materials.resize(2);

	parameters = {
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
}

Settings::Settings(const Configuration &configuration, size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings();
	ParametersReader::configuration(configuration, parameters);
}

Settings::Settings(size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings();
}



