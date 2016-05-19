
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_GENERATOR_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_GENERATOR_H_

#include "../generator.h"
#include "settings.h"
#include "utils.h"

namespace espreso {
namespace input {

template<class TElement>
class CubeGenerator: public UniformGenerator<TElement> {

public:
	static void load(Mesh &mesh, const CubeSettings &settings)
	{
		ESINFO(OVERVIEW) << "Generate cubic mesh";
		ESINFO(DETAILS) << "Cube parameters:\n" << settings;

		CubeGenerator cube(mesh, settings);
		cube.fill();
	}

protected:
	CubeGenerator(Mesh &mesh, const CubeSettings &settings);

	virtual void elementsMaterials(std::vector<Element*> &elements);
	virtual void points(Coordinates &coordinates, size_t &DOFs);
	virtual void boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions);
	virtual void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

	virtual ~CubeGenerator() {};

	const CubeSettings _settings;
	size_t _cluster[3];
};

}
}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_GENERATOR_H_ */
