

#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_

#include "../generator.h"
#include "settings.h"
#include "utils.h"

namespace espreso {
namespace input {

template<class TElement>
class UniformGenerator: public Generator {

public:
	UniformGenerator(Mesh &mesh, const UniformSettings &settings)
	: Generator(mesh, settings), _settings(settings) { };

protected:
	virtual void elementsMesh(std::vector<Element*> &elements);

	virtual void pickElementsInInterval(const std::vector<Element*> &elements, std::vector<Element*> &selection, const Interval &interval);
	virtual void pickNodesInInterval(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const Interval &interval);

	virtual void generateFacesInInterval(std::vector<Element*> &faces, const Interval &interval);
	virtual void generateEdgesInInterval(std::vector<Element*> &edges, const Interval &interval);

	virtual bool partitiate(std::vector<eslocal> &parts);
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints);
	virtual void corners(std::vector<eslocal> &corners);

	TElement _e;
	const UniformSettings _settings;
};

}
}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_ */
