
#include "utils.h"

namespace espreso {
namespace input {

template <class TElement>
eslocal SphereUtils<TElement>::surfaceNodesCount(const SphereSettings &settings)
{
	eslocal cluster[3];
	UniformUtils<TElement>::clusterNodesCount(settings, cluster);

	eslocal surface = 6 * cluster[0] * cluster[1];
	surface -= 4 * (cluster[0] - 2) + 8 * (cluster[1] - 2) + 2 * 8;

	return surface;
}

template <class TElement>
eslocal SphereUtils<TElement>::ringNodesCount(const SphereSettings &settings)
{
	eslocal cluster[3];
	UniformUtils<TElement>::clusterNodesCount(settings, cluster);

	return 4 * cluster[1] - 4;
}

}
}

