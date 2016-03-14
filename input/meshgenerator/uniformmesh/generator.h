

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
	UniformGenerator(const UniformSettings &settings): Generator(settings), _settings(settings) { };

protected:
	virtual void elementsMesh(std::vector<Element*> &elements, std::vector<eslocal> &parts);
	virtual void elementsMaterials(std::vector<Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints);
	virtual void corners(Boundaries &boundaries);

	TElement _e;
	const UniformSettings _settings;
};

}
}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_ */
