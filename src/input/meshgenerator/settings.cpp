
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(Settings &settings)
{
	useMetis = false;

	parameters = {
		{ "USE_METIS" , useMetis   , "Use METIS for mesh partition." },
		{ "MATERIAL1" , material1  , "Parameters of the first material." },
		{ "MATERIAL2" , material2  , "Parameters of the second material." },
		{ "TIME_STEPS", config::solver::TIME_STEPS, "Number of time steps for transient problems."}
	};

	Configuration configuration(Settings::description, options);
}

Settings::Settings(size_t index, size_t size)
: index(index), size(size)
{
	defaultSettings(*this);
}



