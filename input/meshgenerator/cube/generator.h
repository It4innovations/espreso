
#ifndef INPUT_MESHGENERATOR_CUBE_GENERATOR_H_
#define INPUT_MESHGENERATOR_CUBE_GENERATOR_H_

#include "../generator.h"
#include "settings.h"
#include "utils.h"

namespace esinput {

template<class TElement>
class CubeGenerator: public Generator {

public:
	CubeGenerator(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	void fixPoints(std::vector<eslocal> &fixPoints);

	void setForces(mesh::Mesh &mesh, const size_t cluster[]) { };

	void fillGlobalBoundaries(mesh::Boundaries &boundaries, const size_t cluster[]);

	void setCorners(
			mesh::Boundaries &boundaries,
			const size_t number[],
			const bool corners,
			const bool edges,
			const bool surface);

private:
	void fixZeroPlanes(mesh::Mesh &mesh, const size_t cluster[]);
	void fixBottom(mesh::Mesh &mesh, const size_t cluster[]);

	CubeSettings _settings;
	TElement e;
	int _rank;
	int _size;
	size_t _cluster[3];
};

}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_CUBE_GENERATOR_H_ */
