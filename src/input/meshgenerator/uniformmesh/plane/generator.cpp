
#include "generator.h"

namespace espreso {
namespace input {

static void setCluster(size_t cluster[], const PlaneSettings &settings)
{
	if (settings.clusters[0] * settings.clusters[1] != settings.size) {
		ESINFO(espreso::GLOBAL_ERROR) << "The number of clusters(" << settings.clusters[0] * settings.clusters[1]
							<< ") does not accord the number of MPI processes(" << settings.size << ").";
	}

	cluster[2] = 0;
	eslocal index = 0, i = 0;
	for (size_t y = 0; y < settings.clusters[1]; y++) {
		for (size_t x = 0; x < settings.clusters[0]; x++) {
			if (settings.index == index++) {
				cluster[0] = x;
				cluster[1] = y;
				return;
			}
		}
	}
}

template<>
PlaneGenerator<Square4>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Square4>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
PlaneGenerator<Square8>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Square8>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
PlaneGenerator<Triangle3>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Triangle3>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
PlaneGenerator<Triangle6>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Triangle6>(mesh, settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<>
PlaneGenerator<Hexahedron8>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Hexahedron8>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Hexahedron8.";
}

template<>
PlaneGenerator<Hexahedron20>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Hexahedron20>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Hexahedron20.";
}

template<>
PlaneGenerator<Tetrahedron4>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Tetrahedron4>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Tetrahedron4.";
}

template<>
PlaneGenerator<Tetrahedron10>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Tetrahedron10>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Tetrahedron10.";
}

template<>
PlaneGenerator<Prisma6>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Prisma6>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Prisma6.";
}

template<>
PlaneGenerator<Prisma15>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Prisma15>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Prisma15.";
}

template<>
PlaneGenerator<Pyramid5>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Pyramid5>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Pyramid5.";
}

template<>
PlaneGenerator<Pyramid13>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
	: UniformGenerator<Pyramid13>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the element Pyramid13.";
}

}
}
