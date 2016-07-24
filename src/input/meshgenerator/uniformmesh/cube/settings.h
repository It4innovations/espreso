
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

	CubeSettings(const Options &options, size_t index, size_t size);
	CubeSettings(size_t index, size_t size);

	std::vector<Description> description;

	size_t clusters[3];
	double problemLength[3];

	std::map<std::string, std::string> dirichlet;
	std::map<std::string, std::string> forces;
	std::map<std::string, std::string> u;

	std::map<std::string, Interval> regions;

	std::map<std::string, double> heat_sources;

protected:
	void defaultCubeSettings();
};

}
}

#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_SETTINGS_H_ */
