#include "utils.h"

using namespace permoncube;

template <class TElement>
void Utils<TElement>::globalNodesCount(const Settings &settings, esglobal nodes[])
{
	eslocal cluster[3];
	Utils<TElement>::clusterNodesCount(settings, cluster);
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = settings.clusters[i] * (cluster[i] - 1) + 1;
	}
}

template <class TElement>
void Utils<TElement>::clusterNodesCount(const Settings &settings, eslocal nodes[])
{
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] + 1;
	}
}

template <class TElement>
eslocal Utils<TElement>::clusterElementsCount(const Settings &settings)
{
	return TElement::subelements *
	settings.subdomainsInCluster[2] * settings.subdomainsInCluster[1] * settings.subdomainsInCluster[0] *
	settings.elementsInSubdomain[2] * settings.elementsInSubdomain[1] * settings.elementsInSubdomain[0];
}

