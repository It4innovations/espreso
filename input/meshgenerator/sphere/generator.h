
#ifndef INPUT_MESHGENERATOR_SPHERE_GENERATOR_H_
#define INPUT_MESHGENERATOR_SPHERE_GENERATOR_H_

#include "../generator.h"
#include <algorithm>
#include <cmath>
#include "settings.h"
#include "utils.h"

namespace esinput {

template<class TElement>
class SphereGenerator: public Generator {

public:
	SphereGenerator(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	void fixPoints(std::vector<eslocal> &fixPoints);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void corners(mesh::Boundaries &boundaries);
	void clusterBoundaries(mesh::Boundaries &boundaries);

private:
	SphereSettings _settings;
	TElement e;
	int _rank;
	int _size;
};

}

#include "generator.hpp"


#endif /* INPUT_MESHGENERATOR_SPHERE_GENERATOR_H_ */
