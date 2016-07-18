
#include "generator.h"

using namespace espreso::input;

static void setCluster(size_t cluster[], const PlaneSettings &settings)
{
	if (settings.clusters[0] * settings.clusters[1] != settings.size) {
		ESINFO(espreso::GLOBAL_ERROR) << "The number of clusters(" << settings.clusters[0] * settings.clusters[1]
							<< ") does not accord the number of MPI processes(" << settings.size << ").";
	}
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

