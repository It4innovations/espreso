
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_

#include <map>

#include "../settings.h"

namespace espreso {
namespace input {

struct CubeSettings: public UniformSettings {

	enum Faces {
		FRONT,
		REAR,
		LEFT,
		RIGHT,
		TOP,
		BOTTOM
	};

	CubeSettings(const Configuration &configuration, size_t index, size_t size, std::string prefix="");
	CubeSettings(size_t index, size_t size, std::string prefix="");

	std::vector<Description> description;

	size_t clusters[3];
	double problemOrigin[3];
	double problemLength[3];

	std::string projections[3];
	std::string rotations[3];

protected:
	void defaultCubeSettings();
};

}
}

#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_ */
