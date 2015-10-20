
#include "settings.h"

using namespace esinput;

std::vector<Description> createCubeSetting()
{
	std::vector<Description> description(UniformSettings::description);

	std::vector<std::pair<std::string, std::string> > axis = {
			{ "X", "x" },
			{ "Y", "y" },
			{ "Z", "z" }
	};
	std::vector<std::pair<std::string, std::string> > properties = {
			{ "DIRICHLET", "Dirichlet" },
			{ "FORCES", "Force" }
	};
	std::vector<std::pair<std::string, std::string> > cube_faces = {
			{ "FRONT", "front" },
			{ "REAR", "rear" },
			{ "LEFT", "left" },
			{ "RIGHT", "right" },
			{ "TOP", "top" },
			{ "BOTTOM", "bottom" }
	};

	for (size_t i = 0; i < axis.size(); i++) {
		description.push_back({
			INTEGER_PARAMETER, "CLUSTERS_" + axis[i].first, "Number of clusters in clusters in " + axis[i].second + "-axis."
		});
		description.push_back({
			INTEGER_PARAMETER, "LENGTH_" + axis[i].first, "Length of the cube in " + axis[i].second + "-axis."
		});
		for (size_t j = 0; j < properties.size(); j++) {
			for (size_t k = 0; k < cube_faces.size(); k++) {
				description.push_back({
					DOUBLE_PARAMETER,
					properties[j].first + "_" + cube_faces[k].first + "_" + axis[i].first,
					properties[j].second + " on the " + cube_faces[k].second + " face in " + axis[i].second + "-axis."
				});
			}
		}
	}

	return description;
};

std::vector<Description> CubeSettings::description = createCubeSetting();

CubeSettings::CubeSettings(int argc, char** argv): UniformSettings(argc, argv)
{
	Configuration configuration(CubeSettings::description, argc, argv);

	std::vector<std::string> axis = { "X", "Y", "Z" };
	for (size_t i = 0; i < axis.size(); i++) {
		clusters[i] = configuration.value<eslocal>("CLUSTERS_" + axis[i], 1);
		problemLength[i] = configuration.value<double>("LENGTH_" + axis[i], 20);
	}

	std::vector<std::string> properties = { "DIRICHLET", "FORCES" };
	std::vector<std::string> cube_faces = { "FRONT", "REAR", "LEFT", "RIGHT", "TOP", "BOTTOM" };

	fillCondition.resize(cube_faces.size());
	boundaryCondition.resize(cube_faces.size());

	for (size_t f = 0; f < cube_faces.size(); f++) {
		for (size_t p = mesh::DIRICHLET_X; p <= mesh::FORCES_Z; p++) {
			std::string name = properties[p / 3] + "_" + cube_faces[f] + "_" + axis[p % 3];
			fillCondition[f][p] = configuration.isSet(name);
			boundaryCondition[f][p] = configuration.value<double>(name, 0);
		}
	}
}


