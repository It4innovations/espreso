
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_

#include <map>

#include "../settings.h"

namespace espreso {
namespace input {

struct CubeSettings: public UniformSettings {

	enum Faces {
		FRONT,
		REAR,
		LEFT,
		RIGHT,
		TOP,
		BOTTOM
	};

	CubeSettings(const Options &options, size_t index, size_t size);
	CubeSettings(size_t index, size_t size);

	std::vector<Description> description;

	size_t clusters[3];
	double problemLength[3];

	std::vector<double> boundaryCondition;
};

inline std::ostream& operator<<(std::ostream& os, const CubeSettings &s)
{
	os << UniformSettings(s);
	os << "clusters: " << s.clusters[0] << " : " << s.clusters[1] << " : " << s.clusters[2] << "\n";
	os << "cube length: " << s.problemLength[0] << " : " << s.problemLength[1] << " : " << s.problemLength[2] << "\n";

	std::vector<std::string> properties = { "DIRICHLET", "FORCES" };
	std::vector<std::string> cube_faces = { "FRONT", "REAR", "LEFT", "RIGHT", "TOP", "BOTTOM" };
	std::vector<std::string> axis = { "X", "Y", "Z" };

	for (size_t f = 0; f < cube_faces.size(); f++) {
		for (size_t p = DIRICHLET_X; p <= FORCES_Z; p++) {
			std::string name = properties[p / 3] + "_" + cube_faces[f] + "_" + axis[p % 3];
			if (s.boundaryCondition[f * 6 + p] != std::numeric_limits<double>::infinity()) {
				os << name << ": " << s.boundaryCondition[f * 6 + p] << "\n";
			}
		}
	}

	return os;
}

}
}

#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_ */
