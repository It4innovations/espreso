

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

	virtual void partitiate(std::vector<eslocal> &parts);
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints);
	virtual void corners(Boundaries &boundaries);

	TElement _e;
	const UniformSettings _settings;
};

}
}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_ */
