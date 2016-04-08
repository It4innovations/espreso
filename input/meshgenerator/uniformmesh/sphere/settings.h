
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_

#include <map>

#include "../settings.h"

namespace espreso {
namespace input {

struct SphereSettings: public UniformSettings {

	enum Faces {
		INNER,
		OUTER
	};

	SphereSettings(const Options &options, size_t index, size_t size);
	SphereSettings(size_t index, size_t size);

	static std::vector<Description> description;

	size_t layers;
	double innerRadius;
	double outerRadius;

	std::vector<std::map<size_t, double> > boundaryCondition;
	std::vector<std::map<size_t, bool> > fillCondition;
};

inline std::ostream& operator<<(std::ostream& os, const SphereSettings &s)
{
	os << UniformSettings(s);
	os << "layers: " << s.layers << "\n";
	os << "innerRadius: " << s.innerRadius << "\n";
	os << "outerRadius: " << s.outerRadius << "\n";

	std::vector<std::string> properties = { "DIRICHLET", "FORCES" };
	std::vector<std::string> sphere_faces = { "INNER", "OUTER" };
	std::vector<std::string> axis = { "X", "Y", "Z" };

	for (size_t f = 0; f < sphere_faces.size(); f++) {
		for (size_t p = DIRICHLET_X; p <= FORCES_Z; p++) {
			std::string name = properties[p / 3] + "_" + sphere_faces[f] + "_" + axis[p % 3];
			if (s.fillCondition[f].find(p)->second) {
				os << name << ": " << s.boundaryCondition[f].find(p)->second << "\n";
			}
		}
	}

	return os;
}

}
}


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_ */
