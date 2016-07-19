
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

	std::vector<Description> description;

	size_t layers;
	size_t grid;
	double innerRadius;
	double outerRadius;

	std::map<std::string, double> dirichlet;
	std::map<std::string, double> forces;

protected:
	void defaultSphereSettings();
};

}
}


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_SETTINGS_H_ */
