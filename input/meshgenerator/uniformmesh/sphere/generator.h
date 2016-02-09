
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_GENERATOR_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_GENERATOR_H_

#include "../generator.h"
#include "settings.h"
#include "utils.h"

namespace esinput {

template<class TElement>
class SphereGenerator: public UniformGenerator<TElement> {

public:
	SphereGenerator(const SphereSettings &settings);

private:
	void elementsMaterials(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	void points(mesh::Coordinates &coordinates);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Boundaries &boundaries);

	const SphereSettings _settings;
};

}

#include "generator.hpp"


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_GENERATOR_H_ */
