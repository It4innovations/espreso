
#include "settings.h"

#include "esassembler.h"

using namespace espreso::input;

void Settings::defaultSettings()
{
	useMetis = false;

	size_t verbosity = 1;

	parameters = {
		{ prefix + "USE_METIS" , useMetis   , "Use METIS for mesh partition." },
		{ prefix + "MATERIAL1" , material1  , "Parameters of the first material." },
		{ prefix + "MATERIAL2" , material2  , "Parameters of the second material." },
		{ prefix + "TIME_STEPS", config::solver::TIME_STEPS, "Number of time steps for transient problems."},

		{ prefix + "NODES", nodes, "Named sets of nodes.", verbosity },
		{ prefix + "FACES", faces, "Named sets of faces.", verbosity },
		{ prefix + "EDGES", edges, "Named sets of edges.", verbosity },
		{ prefix + "ELEMENTS", elements, "Named sets of elements.", verbosity },

		{ prefix + "DIRICHLET", properties["DIRICHLET"], "Dirichlet boundary conditions.", verbosity },
		{ prefix + "NEUMAN", properties["NEUMAN"], "Neuman boundary conditions.", verbosity },

		{ prefix + "HEAT_SOURCES", properties["HEAT_SOURCES"], "Sources of a heat.", verbosity },
		{ prefix + "TRANSLATION_MOTIONS", properties["TRANSLATION_MOTIONS"], "Translation motion of a region.", verbosity },
		{ prefix + "ACCELERATION", properties["ACCELERATION"], "Acceleration of elements.", verbosity },
		{ prefix + "THICKNESS", properties["THICKNESS"], "Thickness.", verbosity },
		{ prefix + "INITIAL_TEMPERATURE", properties["INITIAL_TEMPERATURE"], "Initial temperature.", verbosity },
		{ prefix + "TEMPERATURE", properties["TEMPERATURE"], "Temperature.", verbosity },
		{ prefix + "OBSTACLE", properties["OBSTACLE"], "Obstacle.", verbosity },
		{ prefix + "NORMAL_DIRECTION", properties["NORMAL_DIRECTION"], "Elements normal direction.", verbosity },

		{ prefix + "INCONSISTENT_STABILIZATION_PARAMETER", AdvectionDiffusion2D::sigma, "Inconsistent stabilization.", verbosity },
		{ prefix + "CONSISTENT_STABILIZATION", AdvectionDiffusion2D::stabilization, "Inconsistent stabilization.", {
				{ "CAU", AdvectionDiffusion2D::STABILIZATION::CAU, "CAU stabilization." },
				{ "SUPG", AdvectionDiffusion2D::STABILIZATION::SUPG, "SUPG stabilization." }
		}, verbosity }
	};
}

Settings::Settings(const ArgsConfiguration &configuration, size_t index, size_t size, std::string prefix)
: index(index), size(size), clusterOffset(0), prefix(prefix)
{
	defaultSettings();
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

Settings::Settings(size_t index, size_t size, std::string prefix)
: index(index), size(size), clusterOffset(0), prefix(prefix)
{
	defaultSettings();
}



