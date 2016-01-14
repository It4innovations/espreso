
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_

#include <map>

#include "../settings.h"

namespace esinput {

struct CubeSettings: public UniformSettings {

	enum Faces {
		FRONT,
		REAR,
		LEFT,
		RIGHT,
		TOP,
		BOTTOM
	};

	CubeSettings(): UniformSettings(0, 1) { };
	CubeSettings(int argc, char** argv, size_t index, size_t size);
	CubeSettings(size_t index, size_t size);

	static std::vector<Description> description;

	size_t clusters[3];
	double problemLength[3];

	std::vector<std::map<size_t, double> > boundaryCondition;
	std::vector<std::map<size_t, bool> > fillCondition;
};

inline std::ostream& operator<<(std::ostream& os, const CubeSettings &s)
{
	os << UniformSettings(s);
	os << "clusters: " << s.clusters[0] << " : " << s.clusters[1] << " : " << s.clusters[2] << "\n";
	os << "cube length: " << s.problemLength[0] << " : " << s.problemLength[1] << " : " << s.problemLength[2] << "\n";
	for (size_t f = 0; f < 6; f++) {
		if (s.fillCondition[f].find(mesh::DIRICHLET_X)->second) {
			os << "DIRICHLET_X: " << s.boundaryCondition[f].find(mesh::DIRICHLET_X)->second << "\n";
		}
		if (s.fillCondition[f].find(mesh::DIRICHLET_Y)->second) {
			os << "DIRICHLET_Y: " << s.boundaryCondition[f].find(mesh::DIRICHLET_Y)->second << "\n";
		}
		if (s.fillCondition[f].find(mesh::DIRICHLET_Z)->second) {
			os << "DIRICHLET_Z: " << s.boundaryCondition[f].find(mesh::DIRICHLET_Z)->second << "\n";
		}
		if (s.fillCondition[f].find(mesh::FORCES_X)->second) {
			os << "FORCES_X: " << s.boundaryCondition[f].find(mesh::FORCES_X)->second << "\n";
		}
		if (s.fillCondition[f].find(mesh::FORCES_Y)->second) {
			os << "FORCES_Y: " << s.boundaryCondition[f].find(mesh::FORCES_Y)->second << "\n";
		}
		if (s.fillCondition[f].find(mesh::FORCES_Z)->second) {
			os << "FORCES_Z: " << s.boundaryCondition[f].find(mesh::FORCES_Z)->second << "\n";
		}
	}
	return os;
}

}

#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_ */
