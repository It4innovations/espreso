
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_

#include "../settings.h"
#include "esmesh.h"

namespace espreso {
namespace input {

struct UniformSettings: public Settings {

	UniformSettings(const ArgsConfiguration &configuration, size_t index, size_t size, std::string prefix="");
	UniformSettings(size_t index, size_t size, std::string prefix="");

	std::vector<Parameter> parameters;

	size_t subdomainsInCluster[3];
	size_t elementsInSubdomain[3];
	size_t materialsLayers[3];

	size_t cornerCount;
	bool corners;
	bool edges;
	bool faces;

protected:
	void defaultUniformSettings();
};

}
}



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_ */
