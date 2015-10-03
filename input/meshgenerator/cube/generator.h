
#ifndef INPUT_MESHGENERATOR_CUBE_GENERATOR_H_
#define INPUT_MESHGENERATOR_CUBE_GENERATOR_H_

#include "settings.h"
#include "../meshgenerator.h"

namespace esinput {

namespace CubeGeneratorOptions {
	enum {
		DIRICHLET_BOTTOM,
		DIRICHLET_ZERO_PLANES
	};
}

template<class TElement>
class CubeGenerator: public MeshGenerator {

public:
	CubeGenerator(const std::string *configFile): e(_settings) {};

	CubeGenerator(CubeSettings &settings): _settings(settings), e(settings)
	{
		_processes = _settings.clusters[0] * _settings.clusters[1] * _settings.clusters[2];
	}

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

#include "../../meshgenerator/cube/generator.hpp"



#endif /* INPUT_MESHGENERATOR_CUBE_GENERATOR_H_ */
