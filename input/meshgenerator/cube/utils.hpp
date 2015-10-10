
#include "utils.h"

using namespace esinput;

template <class TElement>
void CubeUtils<TElement>::globalNodesCount(const CubeSettings &settings, esglobal nodes[])
{
	eslocal cluster[3];
	CubeUtils<TElement>::clusterNodesCount(settings, cluster);
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = settings.clusters[i] * (cluster[i] - 1) + 1;
	}
}

template <class TElement>
void CubeUtils<TElement>::clusterNodesCount(const CubeSettings &settings, eslocal nodes[])
{
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] + 1;
	}
}

template <class TElement>
eslocal CubeUtils<TElement>::clusterElementsCount(const CubeSettings &settings)
{
	return TElement::subelements *
	settings.subdomainsInCluster[2] * settings.subdomainsInCluster[1] * settings.subdomainsInCluster[0] *
	settings.elementsInSubdomain[2] * settings.elementsInSubdomain[1] * settings.elementsInSubdomain[0];
}

