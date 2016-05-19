
#include "utils.h"

namespace espreso {
namespace input {

template <class TElement>
void CubeUtils<TElement>::globalNodesCount(const CubeSettings &settings, esglobal nodes[])
{
	eslocal cluster[3];
	UniformUtils<TElement>::clusterNodesCount(settings, cluster);
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = settings.clusters[i] * (cluster[i] - 1) + 1;
	}
}

}
}
