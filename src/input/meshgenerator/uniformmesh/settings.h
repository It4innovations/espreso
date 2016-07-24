
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_

#include "../settings.h"
#include "esmesh.h"

namespace espreso {
namespace input {

struct UniformSettings: public Settings {

	UniformSettings(const Options &options, size_t index, size_t size);
	UniformSettings(size_t index, size_t size);

	std::vector<Description> description;

	size_t subdomainsInCluster[3];
	size_t elementsInSubdomain[3];
	size_t materialsLayers[3];

	size_t cornerCount;
	bool corners;
	bool edges;
	bool faces;
};

}
}



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_ */
