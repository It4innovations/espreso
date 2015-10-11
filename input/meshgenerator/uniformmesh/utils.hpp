
#include "utils.h"

namespace esinput {


template <class TElement>
void UniformUtils<TElement>::clusterNodesCount(const UniformSettings &settings, eslocal nodes[])
{
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] + 1;
	}
}

template <class TElement>
eslocal UniformUtils<TElement>::clusterElementsCount(const UniformSettings &settings)
{
	return TElement::subelements *
	settings.subdomainsInCluster[2] * settings.subdomainsInCluster[1] * settings.subdomainsInCluster[0] *
	settings.elementsInSubdomain[2] * settings.elementsInSubdomain[1] * settings.elementsInSubdomain[0];
}

}
