
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
			prefix + "CLUSTERS_" + axis[i].first, clusters[i], "Number of clusters in " + axis[i].second + "-axis.", verbosity
		});
	}
	for (size_t i = 0; i < 2; i++) {
		parameters.push_back({
			prefix + "ORIGIN_" + axis[i].first, problemOrigin[i], "Length of the cube in " + axis[i].second + "-axis.", verbosity
		});
	}
	for (size_t i = 0; i < 2; i++) {
		parameters.push_back({
			prefix + "LENGTH_" + axis[i].first, problemLength[i], "Length of the cube in " + axis[i].second + "-axis.", verbosity
		});
	}
	for (size_t i = 0; i < 2; i++) {
		parameters.push_back({
			prefix + "PROJECTION_" + axis[i].first, projections[i], "Projection of " + axis[i].second + "-axis.", verbosity
		});
	}
	parameters.push_back({
		prefix + "ROTATION_" + axis[2].first, rotations[2], "Projection of " + axis[2].second + "-axis.", verbosity
	});

	parameters.push_back({ prefix + "HEAT_SOURCES", properties["HEAT_SOURCES"], "Sources of a heat.", verbosity });
	parameters.push_back({ prefix + "TRANSLATION_MOTIONS", properties["TRANSLATION_MOTIONS"], "Translation motion of a region.", verbosity });
	parameters.push_back({ prefix + "ACCELERATION", properties["ACCELERATION"], "Acceleration of elements.", verbosity });
	parameters.push_back({ prefix + "THICKNESS", properties["THICKNESS"], "Thickness.", verbosity });
	parameters.push_back({ prefix + "INITIAL_TEMPERATURE", properties["INITIAL_TEMPERATURE"], "Initial temperature.", verbosity });
	parameters.push_back({ prefix + "TEMPERATURE", properties["TEMPERATURE"], "Temperature.", verbosity });

	parameters.push_back({ prefix + "INCONSISTENT_STABILIZATION_PARAMETER", AdvectionDiffusion2D::sigma, "Inconsistent stabilization.", verbosity });
	parameters.push_back({ prefix + "CONSISTENT_STABILIZATION", AdvectionDiffusion2D::stabilization, "Inconsistent stabilization.", {
			{ "CAU", AdvectionDiffusion2D::STABILIZATION::CAU, "CAU stabilization." },
			{ "SUPG", AdvectionDiffusion2D::STABILIZATION::SUPG, "SUPG stabilization." }
	}, verbosity });
	parameters.push_back({ prefix + "ELEMENT_BEHAVIOUR", LinearElasticity2D::elementBehaviour, "Element behaviour.", {
			{ "PLAIN_STRAIN", LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRAIN, "Plain strain." },
			{ "PLANE_STRESS", LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS, "Plain stress." },
			{ "PLANE_STRESS_WITH_THICKNESS", LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS, "Plain stress with thickness." },
			{ "AXISYMMETRIC", LinearElasticity2D::ELEMENT_BEHAVIOUR::AXISYMMETRIC, "Axisymmetric." },
	}, verbosity });
}

PlaneSettings::PlaneSettings(const Configuration &configuration, size_t index, size_t size, std::string prefix)
: CubeSettings(index, size, prefix)
{
	parameters.clear();
	defaultPlaneSettings();
	ESINFO(OVERVIEW) << "Load plane setting from file " << configuration.path;
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
	subdomainsInCluster[2] = 1;
	elementsInSubdomain[2] = 1;
}

PlaneSettings::PlaneSettings(size_t index, size_t size, std::string prefix)
: CubeSettings(index, size, prefix)
{
	parameters.clear();
	defaultPlaneSettings();
	parameters.insert(parameters.end(), UniformSettings::parameters.begin(), UniformSettings::parameters.end());
	subdomainsInCluster[2] = 1;
	elementsInSubdomain[2] = 1;
}



