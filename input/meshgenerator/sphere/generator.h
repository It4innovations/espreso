
#ifndef INPUT_MESHGENERATOR_SPHERE_GENERATOR_H_
#define INPUT_MESHGENERATOR_SPHERE_GENERATOR_H_

#include "../generator.h"
#include <algorithm>
#include <cmath>
#include "settings.h"

namespace esinput {

namespace SphereGeneratorOptions {
	enum {
		DIRICHLET_INNER_SURFACE
	};
}

template<class TElement>
class SphereGenerator: public Generator {

public:
	SphereGenerator(int argc, char** argv);

	void fillCluster(int rank, size_t cluster[]);

	void mesh(mesh::Mesh &mesh, const size_t cluster[]);

	void setDirichlet(mesh::Mesh &mesh, const size_t cluster[], size_t dirichlet)
	{
		switch (dirichlet) {
		case SphereGeneratorOptions::DIRICHLET_INNER_SURFACE: {
			fixInnerSurface(mesh, cluster);
			break;
		}
		}
	}

	void setForces(mesh::Mesh &mesh, const size_t cluster[]) { };

	void fillGlobalBoundaries(mesh::Boundaries &boundaries, const size_t cluster[]);

	void setFixPoints(mesh::Mesh &mesh);
	void setCorners(
			mesh::Boundaries &boundaries,
			const size_t number[],
			const bool corners,
			const bool edges,
			const bool surface);

private:
	void fixInnerSurface(mesh::Mesh &mesh, const size_t cluster[]);

	SphereSettings _settings;
	TElement e;
	std::vector<eslocal> _clusterMap;
};

}

#include "../../meshgenerator/sphere/generator.hpp"


#endif /* INPUT_MESHGENERATOR_SPHERE_GENERATOR_H_ */
