
#ifndef SPHERE_H_
#define SPHERE_H_

#include "generator.h"

#include <algorithm>
#include <cmath>

namespace esinput {

namespace SphereGeneratorOptions {
	enum {
		DIRICHLET_INNER_SURFACE
	};
}

struct SphereSettings: public Settings {

	SphereSettings(): Settings()
	{
		layers = 1;
	}

	SphereSettings(Settings &settings): Settings(settings)
	{
		layers = 1;
	}

	size_t layers;
};

inline std::ostream& operator<<(std::ostream& os, const SphereSettings &s)
{
	os << static_cast<Settings>(s);
	os << "layers: " << s.layers << "\n";
	return os;
}

template<class TElement>
class SphereGenerator: public MeshGenerator {

public:
	SphereGenerator(esinput::SphereSettings &settings);

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

#include "sphere.hpp"


#endif /* SPHERE_H_ */
