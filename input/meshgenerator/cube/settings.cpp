
#include "settings.h"

using namespace esinput;

std::vector<Description> CubeSettings::description = {
	{ INTEGER_PARAMETER, "CLUSTERS_X", "Number of clusters in x-axis." },
	{ INTEGER_PARAMETER, "CLUSTERS_Y", "Number of clusters in y-axis." },
	{ INTEGER_PARAMETER, "CLUSTERS_Z", "Number of clusters in z-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_X", "Number of sub-domains in clusters in x-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Y", "Number of sub-domains in clusters in y-axis." },
	{ INTEGER_PARAMETER, "SUBDOMAINS_Z", "Number of sub-domains in clusters in z-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_X", "Number of elements in sub-domains in x-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Y", "Number of elements in sub-domains in y-axis." },
	{ INTEGER_PARAMETER, "ELEMENTS_Z", "Number of elements in sub-domains in z-axis." },
	{ DOUBLE_PARAMETER, "LENGHT_X", "Length of the cube in x-axis." },
	{ DOUBLE_PARAMETER, "LENGHT_Y", "Length of the cube in y-axis." },
	{ DOUBLE_PARAMETER, "LENGHT_Z", "Length of the cube in z-axis." },

	{ DOUBLE_PARAMETER, "DIRICHLET_FRONT_X", "Dirichlet on the front face in x-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_FRONT_Y", "Dirichlet on the front face in y-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_FRONT_Z", "Dirichlet on the front face in z-axis." },

	{ DOUBLE_PARAMETER, "DIRICHLET_REAR_X", "Dirichlet on the rear face in x-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_REAR_Y", "Dirichlet on the rear face in y-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_REAR_Z", "Dirichlet on the rear face in z-axis." },

	{ DOUBLE_PARAMETER, "DIRICHLET_LEFT_X", "Dirichlet on the left face in x-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_LEFT_Y", "Dirichlet on the left face in y-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_LEFT_Z", "Dirichlet on the left face in z-axis." },

	{ DOUBLE_PARAMETER, "DIRICHLET_RIGHT_X", "Dirichlet on the right face in x-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_RIGHT_Y", "Dirichlet on the right face in y-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_RIGHT_Z", "Dirichlet on the right face in z-axis." },

	{ DOUBLE_PARAMETER, "DIRICHLET_TOP_X", "Dirichlet on the top face in x-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_TOP_Y", "Dirichlet on the top face in y-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_TOP_Z", "Dirichlet on the top face in z-axis." },

	{ DOUBLE_PARAMETER, "DIRICHLET_BOTTOM_X", "Dirichlet on the bottom face in x-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_BOTTOM_Y", "Dirichlet on the bottom face in y-axis." },
	{ DOUBLE_PARAMETER, "DIRICHLET_BOTTOM_Z", "Dirichlet on the bottom face in z-axis." },


	{ DOUBLE_PARAMETER, "FORCES_FRONT_X", "Force on the front face in x-axis." },
	{ DOUBLE_PARAMETER, "FORCES_FRONT_Y", "Force on the front face in y-axis." },
	{ DOUBLE_PARAMETER, "FORCES_FRONT_Z", "Force on the front face in z-axis." },

	{ DOUBLE_PARAMETER, "FORCES_REAR_X", "Force on the rear face in x-axis." },
	{ DOUBLE_PARAMETER, "FORCES_REAR_Y", "Force on the rear face in y-axis." },
	{ DOUBLE_PARAMETER, "FORCES_REAR_Z", "Force on the rear face in z-axis." },

	{ DOUBLE_PARAMETER, "FORCES_LEFT_X", "Force on the left face in x-axis." },
	{ DOUBLE_PARAMETER, "FORCES_LEFT_Y", "Force on the left face in y-axis." },
	{ DOUBLE_PARAMETER, "FORCES_LEFT_Z", "Force on the left face in z-axis." },

	{ DOUBLE_PARAMETER, "FORCES_RIGHT_X", "Force on the right face in x-axis." },
	{ DOUBLE_PARAMETER, "FORCES_RIGHT_Y", "Force on the right face in y-axis." },
	{ DOUBLE_PARAMETER, "FORCES_RIGHT_Z", "Force on the right face in z-axis." },

	{ DOUBLE_PARAMETER, "FORCES_TOP_X", "Force on the top face in x-axis." },
	{ DOUBLE_PARAMETER, "FORCES_TOP_Y", "Force on the top face in y-axis." },
	{ DOUBLE_PARAMETER, "FORCES_TOP_Z", "Force on the top face in z-axis." },

	{ DOUBLE_PARAMETER, "FORCES_BOTTOM_X", "Force on the bottom face in x-axis." },
	{ DOUBLE_PARAMETER, "FORCES_BOTTOM_Y", "Force on the bottom face in y-axis." },
	{ DOUBLE_PARAMETER, "FORCES_BOTTOM_Z", "Force on the bottom face in z-axis." }
};

CubeSettings::CubeSettings(int argc, char** argv)
{
	Configuration<CubeSettings> configuration(argc, argv);

	clusters[0] = configuration.value<eslocal>("CLUSTERS_X", 1);
	clusters[1] = configuration.value<eslocal>("CLUSTERS_Y", 1);
	clusters[2] = configuration.value<eslocal>("CLUSTERS_Z", 1);
	subdomainsInCluster[0] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	subdomainsInCluster[1] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	subdomainsInCluster[2] = configuration.value<eslocal>("SUBDOMAINS_X", 2);
	elementsInSubdomain[0] = configuration.value<eslocal>("ELEMENTS_X", 5);
	elementsInSubdomain[1] = configuration.value<eslocal>("ELEMENTS_X", 5);
	elementsInSubdomain[2] = configuration.value<eslocal>("ELEMENTS_X", 5);
	problemLength[0] = configuration.value<double>("LENGTH_X", 20);
	problemLength[1] = configuration.value<double>("LENGTH_Y", 20);
	problemLength[2] = configuration.value<double>("LENGTH_Z", 20);

	std::vector<std::string> properties = { "DIRICHLET", "FORCES" };
	std::vector<std::string> faces = { "FRONT", "REAR", "LEFT", "RIGHT", "TOP", "BOTTOM" };
	std::vector<std::string> axis = { "X", "Y", "Z" };

	fillCondition.resize(faces.size());
	boundaryCondition.resize(faces.size());

	for (size_t f = 0; f < faces.size(); f++) {
		for (size_t p = mesh::DIRICHLET_X; p <= mesh::FORCES_Z; p++) {
			std::string name = properties[p / 3] + "_" + faces[f] + "_" + axis[p % 3];
			fillCondition[f][p] = configuration.isSet(name);
			boundaryCondition[f][p] = configuration.value<double>(name, 0);
		}
	}
}


