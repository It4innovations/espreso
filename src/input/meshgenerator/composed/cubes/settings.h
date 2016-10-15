
#ifndef SRC_INPUT_MESHGENERATOR_COMPOSED_CUBES_SETTINGS_H_
#define SRC_INPUT_MESHGENERATOR_COMPOSED_CUBES_SETTINGS_H_

#include "../../uniformmesh/cube/settings.h"

namespace espreso {
namespace input {

struct CubesSettings {

	CubesSettings(const ArgsConfiguration &configuration, size_t index, size_t size);
	CubesSettings(size_t index, size_t size);

	std::vector<Parameter> parameters;
	CubeSettings cube[2];

protected:
	void defaultCubesSettings();
};

}
}

#endif /* SRC_INPUT_MESHGENERATOR_COMPOSED_CUBES_SETTINGS_H_ */
