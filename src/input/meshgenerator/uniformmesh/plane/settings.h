
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_

#include <map>

#include "../cube/settings.h"

namespace espreso {
namespace input {

struct PlaneSettings: public CubeSettings {

	PlaneSettings(const Configuration &configuration, size_t index, size_t size);
	PlaneSettings(size_t index, size_t size);

protected:
	void defaultPlaneSettings();
};

}
}

#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_SETTINGS_H_ */
