
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
	static void load(Mesh &mesh, const SphereSettings &settings)
	{
		ESINFO(OVERVIEW) << "Generate sphere mesh";
		ParametersReader::printParameters(settings.parameters, config::info::VERBOSE_LEVEL);
		SphereGenerator sphere(mesh, settings);
		sphere.fill();
	}

protected:
	SphereGenerator(Mesh &mesh, const SphereSettings &settings);

	virtual void elementsMaterials(std::vector<Element*> &elements);
	virtual void points(Coordinates &coordinates, size_t &DOFs);
	virtual void boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions);
	virtual void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

	virtual ~SphereGenerator() {};

	const SphereSettings _settings;
	eslocal _cluster[3];
	size_t _side;
};

}
}

#include "generator.hpp"


#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SPHERE_GENERATOR_H_ */
