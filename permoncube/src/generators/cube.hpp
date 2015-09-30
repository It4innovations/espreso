
#ifndef CUBE_HPP_
#define CUBE_HPP_

#include "cube.h"

namespace permoncube {


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

template<class TElement>
void CubeGenerator<TElement>::fillCluster(int rank, size_t cluster[])
{
	eslocal index = 0, i = 0;
	for (size_t z = 0; z < _settings.clusters[2]; z++) {
		for (size_t y = 0; y < _settings.clusters[1]; y++) {
		for (size_t x = 0; x < _settings.clusters[0]; x++) {
				if (rank == index++) {
					cluster[0] = x;
					cluster[1] = y;
					cluster[2] = z;
					return;
				}
			}
		}
	}
}

template <class TElement>
void CubeGenerator<TElement>::mesh(mesh::Mesh &mesh, const size_t cluster[])
{
	for (eslocal i = 0; i < 3; i++) {
		if (_settings.clusters[i] <= cluster[i]) {
			std::cerr << "Permoncube: the number of a cluster is too high\n";
			exit(EXIT_FAILURE);
		}
	}
	eslocal nodes[3];
	Utils<TElement>::clusterNodesCount(_settings, nodes);


	// add coordinates
	mesh::Coordinates &coordinates = mesh.coordinates();
	coordinates.reserve(TElement::clusterNodesCount(_settings));

	esglobal global = 0;
	eslocal local = 0;
	esglobal cs[3], ce[3];
	double step[3];
	for (eslocal i = 0; i < 3; i++) {
		cs[i] = (nodes[i] - 1) * cluster[i];
		ce[i] = (nodes[i] - 1) * (cluster[i] + 1);
	}
	for (eslocal i = 0; i < 3; i++) {
		step[i] = _settings.problemLength[i] / ((nodes[i] - 1) * _settings.clusters[i]);
	}

	esglobal offset[3] = { 0, 0, 0 };
	for (esglobal z = cs[2]; z <= ce[2]; z++) {
		offset[2] = e.offset_z(z);
		for (esglobal y = cs[1]; y <= ce[1]; y++) {
			offset[1] = e.offset_y(y, z);
			for (esglobal x = cs[0]; x <= ce[0]; x++) {
				if (!e.addPoint(x, y, z)) {
					continue;
				}
				offset[0] = e.offset_x(x, y, z);
				global = offset[0] + offset[1] + offset[2];
				coordinates.add(mesh::Point(x * step[0], y * step[1], z * step[2]), local++, global);
			}
		}
	}

	std::vector<eslocal> indices((2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]) * (2 + TElement::subnodes[2]));

	mesh.reserve(Utils<TElement>::clusterElementsCount(_settings));

	eslocal subdomain[3];
	eslocal element[3];

	eslocal subdomainOffset[3];
	eslocal elementOffset[3];

	for (subdomain[2] = 0; subdomain[2] < _settings.subdomainsInCluster[2]; subdomain[2]++) {
		for (subdomain[1] = 0; subdomain[1] < _settings.subdomainsInCluster[1]; subdomain[1]++) {
			for (subdomain[0] = 0; subdomain[0] < _settings.subdomainsInCluster[0]; subdomain[0]++) {
				// for each sub-domain

				for (eslocal i = 0; i < 3; i++) {
					subdomainOffset[i] = subdomain[i] * (_settings.elementsInSubdomain[i] * (1 + TElement::subnodes[i]));
				}
				for (element[2] = 0; element[2] < _settings.elementsInSubdomain[2]; element[2]++) {
					for (element[1] = 0; element[1] < _settings.elementsInSubdomain[1]; element[1]++) {
						for (element[0] = 0; element[0] < _settings.elementsInSubdomain[0]; element[0]++) {
							// for each element

							for (eslocal i = 0; i < 3; i++) {
								elementOffset[i] = subdomainOffset[i] + element[i] * (1 + TElement::subnodes[i]);
							}
							eslocal i = 0;
							for (eslocal z = 0; z < 2 + TElement::subnodes[2]; z++) {
								for (eslocal y = 0; y < 2 + TElement::subnodes[1]; y++) {
									for (eslocal x = 0; x < 2 + TElement::subnodes[0]; x++) {
										// fill node indices

										indices[i++] =
												(elementOffset[2] + z) * nodes[0] * nodes[1] +
												(elementOffset[1] + y) * nodes[0] +
												(elementOffset[0] + x);
									}
								}
							}
							e.addElements(mesh, &indices[0]);
						}
					}
				}
				mesh.endPartition();
			}
		}
	}
}

template <class TElement>
void CubeGenerator<TElement>::fixZeroPlanes(mesh::Mesh &mesh, const size_t cluster[])
{
	mesh::CoordinatesProperty &dirichlet_x = mesh.coordinates().property(mesh::CP::DIRICHLET_X);
	mesh::CoordinatesProperty &dirichlet_y = mesh.coordinates().property(mesh::CP::DIRICHLET_Y);
	mesh::CoordinatesProperty &dirichlet_z = mesh.coordinates().property(mesh::CP::DIRICHLET_Z);

	eslocal nodes[3];
	Utils<TElement>::clusterNodesCount(_settings, nodes);

	if (cluster[0] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal y = 0; y < nodes[1]; y++) {
				if (e.addPoint(0, y, z)) {
					dirichlet_x[e.projectPoint(index)] = 0;
				}
				index += nodes[0];
			}
		}
	}
	if (cluster[1] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (e.addPoint(x, 0, z)) {
					dirichlet_y[e.projectPoint(index)] = 0;
				}
				index++;
			}
			index = (z + 1) * nodes[1] * nodes[0];
		}
	}
	if (cluster[2] == 0) {
		eslocal index = 0;
		for (eslocal y = 0; y < nodes[1]; y++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (e.addPoint(x, y, 0)) {
					dirichlet_z[e.projectPoint(index)] = 0;
				}
				index++;
			}
		}
	}
}

template <class TElement>
void CubeGenerator<TElement>::fixBottom(mesh::Mesh &mesh, const size_t cluster[])
{
	if (cluster[2] > 0) {
		return;
	}
	mesh::CoordinatesProperty &dirichlet_x = mesh.coordinates().property(mesh::CP::DIRICHLET_X);
	mesh::CoordinatesProperty &dirichlet_y = mesh.coordinates().property(mesh::CP::DIRICHLET_Y);
	mesh::CoordinatesProperty &dirichlet_z = mesh.coordinates().property(mesh::CP::DIRICHLET_Z);

	eslocal nodes[3];
	Utils<TElement>::clusterNodesCount(_settings, nodes);

	eslocal index = 0;
	for (eslocal y = 0; y < nodes[1]; y++) {
		for (eslocal x = 0; x < nodes[0]; x++) {
			if (e.addPoint(x, y, 0)) {
				dirichlet_z[e.projectPoint(index)] = 0;
				dirichlet_z[e.projectPoint(index)] = 0;
				dirichlet_y[e.projectPoint(index)] = 0;
				dirichlet_x[e.projectPoint(index)] = 0;
			}
			index++;
		}
	}
}

template <class TElement>
void CubeGenerator<TElement>::fillGlobalBoundaries(mesh::Boundaries &boundaries, const size_t cluster[])
{
	esglobal gNodes[3];
	Utils<TElement>::globalNodesCount(_settings, gNodes);
	eslocal cNodes[3];
	Utils<TElement>::clusterNodesCount(_settings, cNodes);
	boundaries.resize(TElement::clusterNodesCount(_settings));

	bool border[3];
	eslocal cIndex = cluster[0] + cluster[1] * _settings.clusters[0] + cluster[2] * _settings.clusters[0] * _settings.clusters[1];
	esglobal index = 0;

	esglobal cs[3], ce[3];
	for (eslocal i = 0; i < 3; i++) {
		cs[i] = (cNodes[i] - 1) * cluster[i];
		ce[i] = (cNodes[i] - 1) * (cluster[i] + 1);
	}

	for (esglobal z = cs[2]; z <= ce[2]; z++) {
		border[2] = (z == 0 || z == gNodes[2] - 1) ? false : z % ( cNodes[2] - 1) == 0;
		for (esglobal y = cs[1]; y <= ce[1]; y++) {
			border[1] = (y == 0 || y == gNodes[1] - 1) ? false : y % ( cNodes[1] - 1) == 0;
			for (esglobal x = cs[0]; x <= ce[0]; x++) {
				if (!e.addPoint(x, y, z)) {
					continue;
				}
				border[0] = (x == 0 || x == gNodes[0] - 1) ? false : x % ( cNodes[0] - 1) == 0;
				for (int i = 0; i < 8; i++) {
					eslocal tmp = cIndex;
					if (border[0] && (i & 1)) {
						tmp += (x == cs[0]) ? -1 : 1;
					}
					if (border[1] && (i & 2)) {
						tmp += ((y == cs[1]) ? -1 : 1) * _settings.clusters[0];
					}
					if (border[2] && (i & 4)) {
						tmp += ((z == cs[2]) ? -1 : 1) * _settings.clusters[0] * _settings.clusters[1];
					}
					boundaries[index].insert(tmp);
				}
				index++;
			}
		}
	}
}


template <class TElement>
void CubeGenerator<TElement>::setFixPoints(mesh::Mesh &mesh)
{
	std::vector<eslocal> fixPoints;
	fixPoints.reserve(8 * _settings.subdomainsInCluster[0] * _settings.subdomainsInCluster[1] * _settings.subdomainsInCluster[2]);

	eslocal nodes[3];
	eslocal cNodes[3];
	for (int i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
	}
	Utils<TElement>::clusterNodesCount(_settings, cNodes);

	eslocal index;
	eslocal offset[3];
	for (eslocal sz = 0; sz < _settings.subdomainsInCluster[2]; sz++) {
		for (eslocal sy = 0; sy < _settings.subdomainsInCluster[1]; sy++) {
			for (eslocal sx = 0; sx < _settings.subdomainsInCluster[0]; sx++) {
				for (int i = 0; i < 8; i++) {
					offset[0] = (i & 1) ? 1 : 0;
					offset[1] = (i & 2) ? 1 : 0;
					offset[2] = (i & 4) ? 1 : 0;
					index =
							(sz + offset[2]) * nodes[2] * cNodes[0] * cNodes[1] +
							(sy + offset[1]) * nodes[1] * cNodes[0] +
							(sx + offset[0]) * nodes[0];
					fixPoints.push_back(e.projectPoint(index));
				}
			}
		}
	}

	mesh.setFixPoints(fixPoints);
}


template <class TElement>
void CubeGenerator<TElement>::setCorners(
		mesh::Boundaries &boundaries,
		const size_t number[],
		const bool corners,
		const bool edges,
		const bool surface)
{
	eslocal nodes[3];
	eslocal cNodes[3];
	for (int i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
	}
	Utils<TElement>::clusterNodesCount(_settings, cNodes);

	eslocal step[3];
	for (int i = 0; i < 3; i++) {
		step[i] = _settings.elementsInSubdomain[i] / (number[i] + 1);
		step[i] *= TElement::subnodes[i] + 1;
	}
	std::vector<std::vector<size_t> > offsets(3);
	std::vector<size_t> mul(3);

	for (int i = 0; i < 3; i++) {
		for (size_t j = 0; j < _settings.subdomainsInCluster[i]; j++) {
			for (size_t k = 0; k <= number[i] / 2; k++) {
				offsets[i].push_back(j * nodes[i] + k * step[i]);
				offsets[i].push_back(j * nodes[i] + nodes[i] - k * step[i]);
			}
			if (number[i] % 2 == 1) {
				eslocal mid = (_settings.elementsInSubdomain[i] / 2) * (TElement::subnodes[i] + 1);
				offsets[i].push_back(j * nodes[i] + mid);
			}
		}
	}
	mul[0] = 1;
	mul[1] = cNodes[0];
	mul[2] = cNodes[0] * cNodes[1];

	eslocal index;
	for (size_t d = 0; d < 3; d++) {
		for (eslocal i = 1; i < _settings.subdomainsInCluster[d]; i++) {
			for (size_t j = 0; j < offsets[(d + 1) % 3].size(); j++) {
				for (size_t k = 0; k < offsets[(d + 2) % 3].size(); k++) {
					if (!corners
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] == 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] == 0)
					{
						continue;
					}
					if (!edges
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] == 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] != 0)
					{
						continue;
					}
					if (!edges
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] != 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] == 0)
					{
						continue;
					}
					if (!surface
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] != 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] != 0)
					{
						continue;
					}
					index = i * nodes[d] * mul[d];
					index += offsets[(d + 1) % 3][j] * mul[(d + 1) % 3];
					index += offsets[(d + 2) % 3][k] * mul[(d + 2) % 3];
					boundaries.setCorner(e.projectPoint(index));
				}
			}
		}
	}
}

}




#endif /* CUBE_HPP_ */
