
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
		ParametersReader::printParameters(settings.parameters, config::info::VERBOSE_LEVEL);

		CubeGenerator cube(mesh, settings);
		cube.fill();
	}

protected:
	CubeGenerator(Mesh &mesh, const CubeSettings &settings);

	virtual void elementsMaterials(std::vector<Element*> &elements);
	virtual void points(Coordinates &coordinates);
	virtual void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours);

	virtual void settings(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, Coordinates &coordinates);

	virtual ~CubeGenerator() {};

	const CubeSettings _settings;
	size_t _cluster[3];
};

}
}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_GENERATOR_H_ */
