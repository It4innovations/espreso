
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_

#include <map>

#include "../cube/settings.h"

namespace espreso {
namespace input {

struct PlaneSettings: public CubeSettings {

	PlaneSettings(const ArgsConfiguration &configuration, size_t index, size_t size, std::string prefix="");
	PlaneSettings(size_t index, size_t size, std::string prefix="");

protected:
	void defaultPlaneSettings();
};

}
}

#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_ */
