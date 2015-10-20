
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_

#include "../settings.h"

namespace esinput {

struct SphereSettings: public UniformSettings {

	SphereSettings(int argc, char** argv);

	static std::vector<Description> description;

	size_t layers;
	double innerRadius;
	double outerRadius;
};

inline std::ostream& operator<<(std::ostream& os, const SphereSettings &s)
{
	os << "layers: " << s.layers << "\n";
	return os;
}

}


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_ */
