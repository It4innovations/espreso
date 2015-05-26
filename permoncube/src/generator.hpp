#include "generator.h"

using namespace permoncube;


//	###################################################
//	#                                                 #
//	#             A z-coord.                          #
//	#             |                                   #
//	#             |            E3                     #
//	#             |_ _ _ _ _ _ _                      #
//	#            /     E5      /|                     #
//	#           /_ _ _ _ _ _  / |                     #
//	#          |      |      |  |                     #
//	#        E4|      |      |E2|                     #
//	#          |_ _ _ |_ _ _ |  |       y-coord.      #
//	#          |    E1|      |  |------->             #
//	#          |      |      | /                      #
//	#          |_ _ _ |_ _ _ |/                       #
//	#         /                                       #
//	#        /       E0                               #
//	#       /                                         #
//	#      v  x-coord.                                #
//	#                                                 #
//	###################################################

template <class TElement>
void ElementGenerator<TElement>::globalNodesCount(const Settings &settings, size_t nodes[])
{
	size_t cluster[3];
	ElementGenerator<TElement>::clusterNodesCount(settings, cluster);
	for (size_t i = 0; i < 3; i++) {
		nodes[i] = settings.clusters[i] * cluster[i];
		if (settings.clusters[i] > 1) {
			nodes[i]--;
		}
	}
}

template <class TElement>
void ElementGenerator<TElement>::clusterNodesCount(const Settings &settings, size_t nodes[])
{
	for (size_t i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] + 1;
	}
}

template <class TElement>
size_t ElementGenerator<TElement>::clusterElementsCount(const Settings &settings)
{
	return TElement::subelements *
	settings.subdomainsInCluster[2] * settings.subdomainsInCluster[1] * settings.subdomainsInCluster[0] *
	settings.elementsInSubdomain[2] * settings.elementsInSubdomain[1] * settings.elementsInSubdomain[0];
}

template <class TElement>
void ElementGenerator<TElement>::mesh(mesh::Mesh &mesh, const size_t cluster[])
{
	for (int i = 0; i < 3; i++) {
		if (_settings.clusters[i] <= cluster[i]) {
			std::cerr << "Permoncube: the number of a cluster is too high\n";
			exit(EXIT_FAILURE);
		}
	}
	size_t nodes[3];
	globalNodesCount(_settings, nodes);

	TElement::addCoordinates(mesh, _settings, cluster);

	idx_t indices[(2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]) * (2 + TElement::subnodes[2])];

	mesh.reserve(clusterElementsCount(_settings));

	idx_t subdomain[3];
	idx_t element[3];

	idx_t clusterOffset[3];
	idx_t subdomainOffset[3];
	idx_t elementOffset[3];

	for (idx_t i = 0; i < 3; i++) {
		clusterOffset[i] = cluster[i] * (_settings.subdomainsInCluster[i] * (_settings.elementsInSubdomain[i] * (1 + TElement::subnodes[i])));
	}
	for (subdomain[2] = 0; subdomain[2] < _settings.subdomainsInCluster[2]; subdomain[2]++) {
		for (subdomain[1] = 0; subdomain[1] < _settings.subdomainsInCluster[1]; subdomain[1]++) {
			for (subdomain[0] = 0; subdomain[0] < _settings.subdomainsInCluster[0]; subdomain[0]++) {

				// for each sub-domain
				for (idx_t i = 0; i < 3; i++) {
					subdomainOffset[i] = clusterOffset[i] + subdomain[i] * (_settings.elementsInSubdomain[i] * (1 + TElement::subnodes[i]));
				}
				for (element[2] = 0; element[2] < _settings.elementsInSubdomain[2]; element[2]++) {
					for (element[1] = 0; element[1] < _settings.elementsInSubdomain[1]; element[1]++) {
						for (element[0] = 0; element[0] < _settings.elementsInSubdomain[0]; element[0]++) {

							// for each element
							for (idx_t i = 0; i < 3; i++) {
								elementOffset[i] = subdomainOffset[i] + element[i] * (1 + TElement::subnodes[i]);
							}
							idx_t i = 0;
							for (idx_t z = 0; z < 2 + TElement::subnodes[2]; z++) {
								for (idx_t y = 0; y < 2 + TElement::subnodes[1]; y++) {
									for (idx_t x = 0; x < 2 + TElement::subnodes[0]; x++) {
										indices[i++] =
												(elementOffset[2] + z) * nodes[0] * nodes[1] +
												(elementOffset[1] + y) * nodes[0] +
												(elementOffset[0] + x);
									}
								}
							}
							TElement::addElements(mesh, indices);
						}
					}
				}
				mesh.endPartition();
			}
		}
	}
	TElement::clear();
}

