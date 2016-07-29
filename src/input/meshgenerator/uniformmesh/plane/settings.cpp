
#include "settings.h"

#include "esassembler.h"

using namespace espreso::input;

void PlaneSettings::defaultPlaneSettings()
{
	for (size_t i = 0; i < 2; i++) {
		clusters[i] = 1;
		problemLength[i] = 30;
	}

	std::vector<std::pair<std::string, std::string> > axis = {
			{ "X", "x" },
			{ "Y", "y" },
			{ "Z", "z" }
	};

	size_t verbosity = 1;

	for (size_t i = 0; i < 2; i++) {
		parameters.push_back({
			"CLUSTERS_" + axis[i].first, clusters[i], "Number of clusters in " + axis[i].second + "-axis.", verbosity
		});
		parameters.push_back({
			"LENGTH_" + axis[i].first, problemLength[i], "Length of the cube in " + axis[i].second + "-axis.", verbosity
		});
	}

	parameters.push_back({ "NODES", nodes, "Named sets of nodes.", verbosity });
	parameters.push_back({ "FACES", faces, "Named sets of nodes.", verbosity });
	parameters.push_back({ "ELEMENTS", elements, "Named sets of nodes.", verbosity });

	parameters.push_back({ "DIRICHLET", properties["DIRICHLET"], "Dirichlet boundary conditions.", verbosity });
	parameters.push_back({ "HEAT_SOURCES", properties["HEAT_SOURCES"], "Sources of a heat.", verbosity });
	parameters.push_back({ "TRANSLATION_MOTIONS", properties["TRANSLATION_MOTIONS"], "Translation motion of a region.", verbosity });

	parameters.push_back({ "INCONSISTENT_STABILIZATION_PARAMETER", AdvectionDiffusion2D::sigma, "Inconsistent stabilization.", verbosity });
	parameters.push_back({ "CONSISTENT_STABILIZATION", AdvectionDiffusion2D::stabilization, "Inconsistent stabilization.", {
			{ "CAU", AdvectionDiffusion2D::STABILIZATION::CAU, "CAU stabilization." },
			{ "SUPG", AdvectionDiffusion2D::STABILIZATION::SUPG, "SUPG stabilization." }
	}, verbosity });
}

PlaneSettings::PlaneSettings(const Configuration &configuration, size_t index, size_t size)
: CubeSettings(index, size)
{
	parameters.clear();
	defaultPlaneSettings();
	ESINFO(OVERVIEW) << "Load plane setting from file " << configuration.path;
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
	subdomainsInCluster[2] = 1;
	elementsInSubdomain[2] = 1;
}

PlaneSettings::PlaneSettings(size_t index, size_t size)
: CubeSettings(index, size)
{
	parameters.clear();
	defaultPlaneSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	subdomainsInCluster[2] = 1;
	elementsInSubdomain[2] = 1;
}


