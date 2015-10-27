

#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_

#include "../generator.h"
#include "settings.h"
#include "utils.h"

namespace esinput {

template<class TElement>
class UniformGenerator: public Generator {

public:
	UniformGenerator(int argc, char** argv, int rank, int size);

	virtual void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints);
	virtual void corners(mesh::Boundaries &boundaries);

protected:
	TElement _e;
	UniformSettings _settings;
};

}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_ */
