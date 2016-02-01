

#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_

#include "../generator.h"
#include "settings.h"
#include "utils.h"

namespace esinput {

template<class TElement>
class UniformGenerator: public Generator {

protected:
	UniformGenerator(int argc, char** argv, size_t index, size_t size)
	: Generator(argc, argv, index, size), _settings(argc, argv, index, size) { };
	UniformGenerator(const UniformSettings &settings)
	: Generator(settings), _settings(settings) { };

	virtual void elementsMesh(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	virtual void elementsMaterials(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints);
	virtual void corners(mesh::Boundaries &boundaries);

	TElement _e;
	const UniformSettings _settings;
};

}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_GENERATOR_H_ */
