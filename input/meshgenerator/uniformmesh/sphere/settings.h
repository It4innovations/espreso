
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_

#include "../settings.h"

namespace esinput {

struct SphereSettings: public UniformSettings {

	SphereSettings(const Options &options, size_t index, size_t size);
	SphereSettings(size_t index, size_t size);

	static std::vector<Description> description;

	size_t layers;
	double innerRadius;
	double outerRadius;
};

inline std::ostream& operator<<(std::ostream& os, const SphereSettings &s)
{
	os << UniformSettings(s);
	os << "layers: " << s.layers << "\n";
	os << "innerRadius: " << s.innerRadius << "\n";
	os << "outerRadius: " << s.outerRadius << "\n";
	return os;
}

}


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_ */
