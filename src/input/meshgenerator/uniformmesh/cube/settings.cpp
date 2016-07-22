
#include "settings.h"

using namespace espreso::input;

static void defaultSettings(CubeSettings &settings)
{
	for (size_t i = 0; i < 3; i++) {
		settings.clusters[i] = 1;
		settings.problemLength[i] = 30;
	}

	std::vector<std::pair<std::string, std::string> > axis = {
			{ "X", "x" },
			{ "Y", "y" },
			{ "Z", "z" }
	};

	size_t verbosity = 1;

	for (size_t i = 0; i < axis.size(); i++) {
		parameters.push_back({
			"CLUSTERS_" + axis[i].first, clusters[i], "Number of clusters in " + axis[i].second + "-axis.", verbosity
		});
		parameters.push_back({
			"LENGTH_" + axis[i].first, problemLength[i], "Length of the cube in " + axis[i].second + "-axis.", verbosity
		});
	}

	parameters.push_back({ "DIRICHLET", dirichlet, "Dirichlet boundary conditions.", verbosity });
	parameters.push_back({ "FORCES", forces, "Boundary forces.", verbosity });

	parameters.push_back({ "REGIONS", regions, "Named sets of nodes.", verbosity });
}

	Configuration configuration(CubeSettings::description, options);
}

CubeSettings::CubeSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	defaultSettings(*this);
}


