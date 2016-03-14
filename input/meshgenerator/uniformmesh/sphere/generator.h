
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_GENERATOR_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_GENERATOR_H_

#include "../generator.h"
#include "settings.h"
#include "utils.h"

namespace espreso {
namespace input {

template<class TElement>
class SphereGenerator: public UniformGenerator<TElement> {

public:
	SphereGenerator(const SphereSettings &settings);

private:
	void elementsMaterials(std::vector<Element*> &elements, std::vector<eslocal> &parts);
	void points(Coordinates &coordinates);
	void boundaryConditions(Coordinates &coordinates);
	void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

	const SphereSettings _settings;
};

}
}

#include "generator.hpp"


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_GENERATOR_H_ */
