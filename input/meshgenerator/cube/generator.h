
#ifndef INPUT_MESHGENERATOR_CUBE_GENERATOR_H_
#define INPUT_MESHGENERATOR_CUBE_GENERATOR_H_

#include "../generator.h"
#include "settings.h"

namespace esinput {

namespace CubeGeneratorOptions {
	enum {
		DIRICHLET_BOTTOM,
		DIRICHLET_ZERO_PLANES
	};
}

template<class TElement>
class CubeGenerator: public Generator {

public:
	CubeGenerator(int argc, char** argv): _settings(argc, argv), e(_settings) {};

	void fillCluster(int rank, size_t cluster[]);

	void mesh(mesh::Mesh &mesh, const size_t cluster[]);

	void setDirichlet(mesh::Mesh &mesh, const size_t cluster[], size_t dirichlet)
	{
		switch (dirichlet) {
		case CubeGeneratorOptions::DIRICHLET_BOTTOM: {
			fixBottom(mesh, cluster);
			break;
		}
		case CubeGeneratorOptions::DIRICHLET_ZERO_PLANES: {
			fixZeroPlanes(mesh, cluster);
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
	void fixZeroPlanes(mesh::Mesh &mesh, const size_t cluster[]);
	void fixBottom(mesh::Mesh &mesh, const size_t cluster[]);

	CubeSettings _settings;
	TElement e;
};

}

#include "generator.hpp"



#endif /* INPUT_MESHGENERATOR_CUBE_GENERATOR_H_ */
