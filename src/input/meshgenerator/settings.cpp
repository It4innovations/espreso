
#include "settings.h"

using namespace espreso::input;

void Settings::defaultSettings()
{
	useMetis    = false;

	materials.resize(2);

	parameters = {
		{ "USE_METIS"   , useMetis   , "Use METIS for mesh partition." },

		{ "MAT1_DENSITY", materials[0].density     , "Density of the first material." },
		{ "MAT1_YOUNG"  , materials[0].youngModulus, "Young's modulus of the first material." },
		{ "MAT1_POISSON", materials[0].poissonRatio, "Poisson's ratio of the first material." },
		{ "MAT2_DENSITY", materials[1].density     , "Density of the first material." },
		{ "MAT2_YOUNG"  , materials[1].youngModulus, "Young's modulus of the first material." },
		{ "MAT2_POISSON", materials[1].poissonRatio, "Poisson's ratio of the first material." },

		{ "TIME_STEPS", config::solver::TIME_STEPS, "Number of time steps for transient problems."}
	};
}

Settings::Settings(const Configuration &configuration, size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings();
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

Settings::Settings(size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings();
}



