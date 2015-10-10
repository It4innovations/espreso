#include "utils.h"

using namespace esinput;

template <class TElement>
void SphereUtils<TElement>::clusterNodesCount(const SphereSettings &settings, eslocal nodes[])
{
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] + 1;
	}
}

template <class TElement>
eslocal SphereUtils<TElement>::clusterElementsCount(const SphereSettings &settings)
{
	return TElement::subelements *
	settings.subdomainsInCluster[2] * settings.subdomainsInCluster[1] * settings.subdomainsInCluster[0] *
	settings.elementsInSubdomain[2] * settings.elementsInSubdomain[1] * settings.elementsInSubdomain[0];
}

