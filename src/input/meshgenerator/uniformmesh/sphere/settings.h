
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

	SphereSettings(const ArgsConfiguration &configuration, size_t index, size_t size, std::string prefix="");
	SphereSettings(size_t index, size_t size, std::string prefix="");

	std::vector<Parameter> parameters;

	size_t layers;
	size_t grid;
	double innerRadius;
	double outerRadius;

protected:
	void defaultSphereSettings();
};

}
}


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_ */
