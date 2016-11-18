
#include "generator.h"

namespace espreso {
namespace input {

static void checkSettings(size_t cluster[], size_t &side, const SphereSettings &settings)
{
	if (settings.grid * settings.grid * settings.layers * 6 != settings.size) {
		ESINFO(espreso::GLOBAL_ERROR) << "The number of clusters(" << settings.grid * settings.grid * settings.layers * 6
				<< ") does not accord the number of MPI processes(" << settings.size << ").";
	}

	if (settings.subdomainsInCluster[1] * settings.elementsInSubdomain[1] != settings.subdomainsInCluster[0] * settings.elementsInSubdomain[0]) {
		ESINFO(espreso::GLOBAL_ERROR) << "The number of elements in x-axis does not accord the number of elements in y-axis";
	}

	side = settings.index / (settings.size / 6);

	size_t index = 0, sideIndex = settings.index % (settings.size / 6);

	for (size_t z = 0; z < settings.layers; z++) {
		for (size_t y = 0; y < settings.grid; y++) {
			for (size_t x = 0; x < settings.grid; x++) {
				if (sideIndex == index++) {
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
SphereGenerator<Hexahedron8>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Hexahedron8>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}

template<>
SphereGenerator<Hexahedron20>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Hexahedron20>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}

template<>
SphereGenerator<Tetrahedron4>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Tetrahedron4>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}

template<>
SphereGenerator<Tetrahedron10>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Tetrahedron10>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}

template<>
SphereGenerator<Prisma6>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Prisma6>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}

template<>
SphereGenerator<Prisma15>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Prisma15>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}

template<>
SphereGenerator<Pyramid5>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Pyramid5>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}

template<>
SphereGenerator<Pyramid13>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<Pyramid13>(mesh, settings), _settings(settings)
{
	checkSettings(_cluster, _side, _settings);
}


}
}
