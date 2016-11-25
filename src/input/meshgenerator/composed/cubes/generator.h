
#ifndef SRC_INPUT_MESHGENERATOR_COMPOSED_CUBES_GENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_COMPOSED_CUBES_GENERATOR_H_

#include "../../uniformmesh/cube/generator.h"
#include "settings.h"

namespace espreso {
namespace input {

template<class TElement>
class CubesGenerator: public Loader {

public:
	static void load(Mesh &mesh, CubesSettings &settings)
	{
		ParametersReader::printParameters(settings.parameters, config::info::VERBOSE_LEVEL);

		CubesGenerator cubes(mesh, settings);
		ESINFO(OVERVIEW) << "Generate cubes connected by Mortar interface";
		cubes.fill();

	}

protected:
	CubesGenerator(Mesh &mesh, CubesSettings &settings);

	virtual void points(Coordinates &coordinates);
	virtual void elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	virtual void materials(std::vector<Material> &materials);
	virtual void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours);
	virtual void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);

	virtual bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners);

	virtual ~CubesGenerator() {};

	CubesSettings _settings;
	size_t _cluster[3];

	size_t _cubeIndex;
	Loader* _loader;
	TElement _e;
};

}
}

#include "generator.hpp"



#endif /* SRC_INPUT_MESHGENERATOR_COMPOSED_CUBES_GENERATOR_H_ */
