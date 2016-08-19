
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(Settings &settings)
{
	useMetis = false;

	size_t verbosity = 1;

	parameters = {
		{ prefix + "USE_METIS" , useMetis   , "Use METIS for mesh partition." },
		{ prefix + "MATERIAL1" , material1  , "Parameters of the first material." },
		{ prefix + "MATERIAL2" , material2  , "Parameters of the second material." },
		{ prefix + "TIME_STEPS", config::solver::TIME_STEPS, "Number of time steps for transient problems."},

		{ prefix + "NODES", nodes, "Named sets of nodes.", verbosity },
		{ prefix + "FACES", faces, "Named sets of nodes.", verbosity },
		{ prefix + "ELEMENTS", elements, "Named sets of nodes.", verbosity },

		{ prefix + "DIRICHLET", properties["DIRICHLET"], "Dirichlet boundary conditions.", verbosity }
	};
}

Settings::Settings(const Configuration &configuration, size_t index, size_t size, std::string prefix)
: index(index), size(size), clusterOffset(0), prefix(prefix)
{
	defaultSettings();
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

Settings::Settings(size_t index, size_t size, std::string prefix)
: index(index), size(size), clusterOffset(0), prefix(prefix)
{
	defaultSettings(*this);
}



