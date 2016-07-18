
#include "generator.h"

using namespace espreso::input;

static void setCluster(size_t cluster[], const CubeSettings &settings)
{
	if (settings.clusters[0] * settings.clusters[1] * settings.clusters[2] != settings.size) {
		ESINFO(espreso::GLOBAL_ERROR) << "The number of clusters(" << settings.clusters[0] * settings.clusters[1] * settings.clusters[2]
							<< ") does not accord the number of MPI processes(" << settings.size << ").";
	}
	eslocal index = 0, i = 0;
	for (size_t z = 0; z < settings.clusters[2]; z++) {
		for (size_t y = 0; y < settings.clusters[1]; y++) {
			for (size_t x = 0; x < settings.clusters[0]; x++) {
				if (settings.index == index++) {
					cluster[0] = x;
					cluster[1] = y;
					cluster[2] = z;
					return;
				}
			}
		}
	}
}

template<>
CubeGenerator<Hexahedron8>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Hexahedron8>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
CubeGenerator<Hexahedron20>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Hexahedron20>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
CubeGenerator<Tetrahedron4>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Tetrahedron4>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
CubeGenerator<Tetrahedron10>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Tetrahedron10>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
CubeGenerator<Prisma6>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Prisma6>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
CubeGenerator<Prisma15>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Prisma15>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
CubeGenerator<Pyramid5>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Pyramid5>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
CubeGenerator<Pyramid13>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
	: UniformGenerator<Pyramid13>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}



