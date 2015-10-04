
#ifndef CUBE_H_
#define CUBE_H_

#include "generator.h"

namespace esinput {

namespace CubeGeneratorOptions {
	enum {
		DIRICHLET_BOTTOM,
		DIRICHLET_ZERO_PLANES
	};
}

struct CubeSettings: public Settings {

	CubeSettings(): Settings()
	{
		problemLength[0] = problemLength[1] = problemLength[2] = 30;
	}

	CubeSettings(Settings &settings): Settings(settings)
	{
		problemLength[0] = problemLength[1] = problemLength[2] = 30;
	};

	double problemLength[3];
};

inline std::ostream& operator<<(std::ostream& os, const CubeSettings &s)
{
	os << static_cast<Settings>(s);
	os << "clusterLength: " << s.problemLength[0] << " : " << s.problemLength[1] << " : " << s.problemLength[2] << "\n";
	return os;
}

template<class TElement>
class CubeGenerator: public Generator {

public:
	CubeGenerator(esinput::CubeSettings &settings): _settings(settings), e(settings)
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

#include "cube.hpp"



#endif /* CUBE_H_ */
