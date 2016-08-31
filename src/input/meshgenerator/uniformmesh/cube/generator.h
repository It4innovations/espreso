
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

	CubeGenerator(Mesh &mesh, const CubeSettings &settings);
protected:
	virtual void elementsMaterials(std::vector<Element*> &elements);

	virtual void pickElementsInInterval(const std::vector<Element*> &elements, std::vector<Element*> &selection, const Interval &interval);
	virtual void pickNodesInInterval(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const Interval &interval);
	virtual void generateFacesInInterval(std::vector<Element*> &faces, const Interval &interval);
	virtual void generateEdgesInInterval(std::vector<Element*> &edges, const Interval &interval);

	virtual void points(Coordinates &coordinates);
	virtual void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours);

	virtual ~CubeGenerator() {};

	const CubeSettings _settings;
	size_t _cluster[3];
};

}
}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_CUBE_GENERATOR_H_ */
