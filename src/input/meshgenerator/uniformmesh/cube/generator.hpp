
#include "generator.h"

namespace espreso {
namespace input {


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

// Intel compilator error: do not use generic contructor

//template<class TElement>
//CubeGenerator<TElement>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
//: UniformGenerator<TElement>(mesh, settings), _settings(settings)
//{
//	ESINFO(GLOBAL_ERROR) << "Cube generator does not support the selected element type.";
//}

template<class TElement>
void CubeGenerator<TElement>::elementsMaterials(std::vector<Element*> &elements)
{
	esglobal cubeElements[3], partSize[3], cOffset[3], offset[3];
	eslocal subdomain[3], element[3], material, counter;

	for (size_t i = 0; i < 3; i++) {
		cubeElements[i] = _settings.clusters[i] * _settings.subdomainsInCluster[i] * _settings.elementsInSubdomain[i];
		cOffset[i] = _cluster[i] * _settings.subdomainsInCluster[i] * _settings.elementsInSubdomain[i];
		partSize[i] = std::ceil(cubeElements[i] / (double)_settings.materialsLayers[i]);
	}

	counter = 0;
	for (subdomain[2] = 0; subdomain[2] < _settings.subdomainsInCluster[2]; subdomain[2]++) {
			for (subdomain[1] = 0; subdomain[1] < _settings.subdomainsInCluster[1]; subdomain[1]++) {
				for (subdomain[0] = 0; subdomain[0] < _settings.subdomainsInCluster[0]; subdomain[0]++) {

					for (element[2] = 0; element[2] < _settings.elementsInSubdomain[2]; element[2]++) {
						for (element[1] = 0; element[1] < _settings.elementsInSubdomain[1]; element[1]++) {
							for (element[0] = 0; element[0] < _settings.elementsInSubdomain[0]; element[0]++) {

								material = 0;
								for (eslocal i = 0; i < 3; i++) {
									offset[i] = cOffset[i] + subdomain[i] * _settings.elementsInSubdomain[i] + element[i];
									if (offset[i] / partSize[i] % 2 == 1) {
										material = (material + 1) % 2;
									}
								}
								for (size_t e = 0; e < TElement::subelements; e++) {
									elements[counter++]->setParam(Element::MATERIAL, material);
								}
							}
						}
					}

				}
			}
	}

}


template<class TElement>
void CubeGenerator<TElement>::points(Coordinates &coordinates)
{
	eslocal cNodes[3];
	esglobal gNodes[3];

	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	CubeUtils<TElement>::globalNodesCount(_settings, gNodes);

	coordinates.clear();
	coordinates.reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	esglobal cs[3], ce[3];
	double step[3];
	for (eslocal i = 0; i < 3; i++) {
		cs[i] = (cNodes[i] - 1) * _cluster[i];
		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
	}
	for (eslocal i = 0; i < 3; i++) {
		step[i] = _settings.problemLength[i] / ((cNodes[i] - 1) * _settings.clusters[i]);
	}

	for (esglobal z = cs[2]; z <= ce[2]; z++) {
		for (esglobal y = cs[1]; y <= ce[1]; y++) {
			for (esglobal x = cs[0]; x <= ce[0]; x++) {
				coordinates.add(
					Point(_settings.problemOrigin[0] + x * step[0], _settings.problemOrigin[1] + y * step[1], _settings.problemOrigin[2] + z * step[2]),
					(z - cs[2]) * cNodes[0] * cNodes[1] + (y - cs[1]) * cNodes[0] + (x - cs[0]),
					z * gNodes[0] * gNodes[1] + y * gNodes[0] + x
				);
			}
		}
	}
}


template<class TElement>
static void goThroughElements(
		const CubeSettings &settings,
		const Interval &interval,
		size_t cluster[],
		std::function<void(std::vector<eslocal> &indices, CubeEdges edge)> addEdge,
		std::function<void(std::vector<eslocal> &indices, CubeFaces face)> addFace,
		bool restrictNodes)
{
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(settings, cNodes);
	std::vector<eslocal> indices((2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]) * (2 + TElement::subnodes[2]));

	size_t start[3], end[3];
	CubeUtils<TElement>::computeInterval(settings, cluster, interval, start, end);

	CubeEdges edge = CubeUtils<TElement>::cubeEdge(settings, cluster, interval);
	CubeFaces face = CubeUtils<TElement>::cubeFace(settings, cluster, interval);

	if (edge == CubeEdges::NONE && face == CubeFaces::NONE) {
		return;
	}

	size_t minOffset[3], maxOffset[3];
	for (size_t i = 0; i < 3; i++) {
		double min = settings.clusters[i] * (cNodes[i] - 1) * interval.start[i] / settings.problemLength[i];
		double max = settings.clusters[i] * (cNodes[i] - 1) * interval.end[i] / settings.problemLength[i];
		if (min == std::round(min)) {
			minOffset[i] = std::round(min) + (interval.excludeStart[i] ? 1 : 0);
		} else {
			minOffset[i] = std::ceil(min);
		}
		if (max == std::round(max)) {
			maxOffset[i] = std::round(max) - (interval.excludeEnd[i] ? 1 : 0);
		} else {
			maxOffset[i] = std::floor(max);
		}
		if (minOffset[i] > maxOffset[i]) {
			return;
		}

		size_t cStart = cluster[i] * (cNodes[i] - 1);
		size_t cEnd = (cluster[i] + 1) * (cNodes[i] - 1);

		if (maxOffset[i] < cStart || cEnd < minOffset[i]) {
			return;
		}
		minOffset[i] = minOffset[i] < cStart ? 0 : minOffset[i] - cStart;
		maxOffset[i] = maxOffset[i] > cEnd ? cNodes[i] - 1 : maxOffset[i] - cStart;
	}

	for (size_t ex = start[0]; ex < end[0]; ex++) {
		for (size_t ey = start[1]; ey < end[1]; ey++) {
			for (size_t ez = start[2]; ez < end[2]; ez++) {

				for (eslocal z = 0, i = 0; z <= 1 + TElement::subnodes[2]; z++) {
					for (eslocal y = 0; y <= 1 + TElement::subnodes[1]; y++) {
						for (eslocal x = 0; x <= 1 + TElement::subnodes[0]; x++, i++) {


							size_t offsetX = ex * (1 + TElement::subnodes[0]) + x;
							size_t offsetY = ey * (1 + TElement::subnodes[1]) + y;
							size_t offsetZ = ez * (1 + TElement::subnodes[2]) + z;

							if (restrictNodes) {
								offsetX = offsetX < minOffset[0] ? minOffset[0] : maxOffset[0] < offsetX ? maxOffset[0] : offsetX;
								offsetY = offsetY < minOffset[1] ? minOffset[1] : maxOffset[1] < offsetY ? maxOffset[1] : offsetY;
								offsetZ = offsetZ < minOffset[2] ? minOffset[2] : maxOffset[2] < offsetZ ? maxOffset[2] : offsetZ;
							}

							indices[i] = offsetZ * cNodes[1] * cNodes[0] + offsetY * cNodes[0] + offsetX;
						}
					}
				}
				if (edge != CubeEdges::NONE) {
					addEdge(indices, edge);
				}
				if (face != CubeFaces::NONE) {
					addFace(indices, face);
				}
			}
		}
	}
}

template<class TElement>
void CubeGenerator<TElement>::pickElementsInInterval(const std::vector<Element*> &elements, std::vector<Element*> &selection, const Interval &interval)
{
	ESINFO(GLOBAL_ERROR) << "Implement pick elements in interval";
}

template<class TElement>
void CubeGenerator<TElement>::pickNodesInInterval(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const Interval &interval)
{
	goThroughElements<TElement>(
			_settings, interval, _cluster,
			[ & ] (std::vector<eslocal> &indices, CubeEdges edge) {
				this->_e.pickNodes(nodes, selection, indices.data(), edge);
			},
			[ & ] (std::vector<eslocal> &indices, CubeFaces face) {
				this->_e.pickNodes(nodes, selection, indices.data(), face);
			},
			true
	);

	std::sort(selection.begin(), selection.end());
	Esutils::removeDuplicity(selection);
}

template<class TElement>
void CubeGenerator<TElement>::generateFacesInInterval(std::vector<Element*> &faces, const Interval &interval)
{
	goThroughElements<TElement>(
			_settings, interval, _cluster,
			[ & ] (std::vector<eslocal> &indices, CubeEdges edge) {
				ESINFO(GLOBAL_ERROR) << "Invalid interval";
			},
			[ & ] (std::vector<eslocal> &indices, CubeFaces face) {
				this->_e.addFaces(faces, &indices[0], face);
			},
			false
	);
}

template<class TElement>
void CubeGenerator<TElement>::generateEdgesInInterval(std::vector<Element*> &edges, const Interval &interval)
{
	goThroughElements<TElement>(
			_settings, interval, _cluster,
			[ & ] (std::vector<eslocal> &indices, CubeEdges edge) {
				this->_e.addEdges(edges, &indices[0], edge);
			},
			[ & ] (std::vector<eslocal> &indices, CubeFaces face) {
				ESINFO(GLOBAL_ERROR) << "Implement add edges on face";
			},
			false
	);
}


template <class TElement>
void CubeGenerator<TElement>::clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours)
{
	esglobal gNodes[3];
	CubeUtils<TElement>::globalNodesCount(_settings, gNodes);
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	bool border[3];
	eslocal cIndex = _cluster[0] + _cluster[1] * _settings.clusters[0] + _cluster[2] * _settings.clusters[0] * _settings.clusters[1];
	esglobal index = 0;

	esglobal cs[3], ce[3];
	for (eslocal i = 0; i < 3; i++) {
		cs[i] = (cNodes[i] - 1) * _cluster[i];
		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
	}

	// TODO: optimize this
	std::set<int> neighs;

	for (esglobal z = cs[2]; z <= ce[2]; z++) {
		border[2] = (z == 0 || z == gNodes[2] - 1) ? false : z % ( cNodes[2] - 1) == 0;
		for (esglobal y = cs[1]; y <= ce[1]; y++) {
			border[1] = (y == 0 || y == gNodes[1] - 1) ? false : y % ( cNodes[1] - 1) == 0;
			for (esglobal x = cs[0]; x <= ce[0]; x++) {
				border[0] = (x == 0 || x == gNodes[0] - 1) ? false : x % ( cNodes[0] - 1) == 0;
				for (int i = 0; i < 8; i++) {
					eslocal tmp = cIndex + _settings.clusterOffset;
					if (border[0] && (i & 1)) {
						tmp += (x == cs[0]) ? -1 : 1;
					}
					if (border[1] && (i & 2)) {
						tmp += ((y == cs[1]) ? -1 : 1) * _settings.clusters[0];
					}
					if (border[2] && (i & 4)) {
						tmp += ((z == cs[2]) ? -1 : 1) * _settings.clusters[0] * _settings.clusters[1];
					}
					nodes[index]->clusters().push_back(tmp);
					neighs.insert(tmp);
				}
				std::sort(nodes[index]->clusters().begin(), nodes[index]->clusters().end());
				Esutils::removeDuplicity(nodes[index]->clusters());
				index++;
			}
		}
	}

	neighs.erase(config::env::MPIrank);
	neighbours.insert(neighbours.end(), neighs.begin(), neighs.end());
}

}
}

