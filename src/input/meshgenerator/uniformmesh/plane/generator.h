
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_GENERATOR_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_GENERATOR_H_

#include "../cube/generator.h"
#include "settings.h"

namespace espreso {
namespace input {

template<class TElement>
class PlaneGenerator: public UniformGenerator<TElement> {

public:
	static void load(Mesh &mesh, const PlaneSettings &settings)
	{
		ESINFO(OVERVIEW) << "Generate plane mesh";
		ParametersReader::printParameters(settings.parameters, config::info::VERBOSE_LEVEL);

		PlaneGenerator plane(mesh, settings);
		plane.fill();
	}

protected:
	PlaneGenerator(Mesh &mesh, const PlaneSettings &settings);

	virtual void elementsMesh(std::vector<Element*> &elements);
	virtual void elementsMaterials(std::vector<Element*> &elements);

	virtual void generateFacesInInterval(std::vector<Element*> &faces, const Interval &interval);
	virtual void generateEdgesInInterval(std::vector<Element*> &edges, const Interval &interval);
	virtual void generateNodesInInterval(std::vector<Element*> &nodes, const Interval &interval);

	virtual void points(Coordinates &coordinates);
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints);
	virtual void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours);
	virtual void corners(std::vector<eslocal> &corners);

	virtual void settings(
			std::vector<Evaluator*> &evaluators,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);

	virtual ~PlaneGenerator() {};

	const PlaneSettings _settings;
	size_t _cluster[3];
};

}
}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_PLANE_GENERATOR_H_ */
