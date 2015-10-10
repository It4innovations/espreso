
#include "generator.h"

namespace esinput {


//	###################################################
//	#                                                 #
//	#             A z-coord.                          #
//	#             |                                   #
//	#             |                                   #
//	#             |                                   #
//	#             |                                   #
//	#             |                                   #
//	#             |                                   #
//	#             |                                   #
//	#             |                y-coord.           #
//	#             /--------------------->             #
//	#            /                                    #
//	#           /                                     #
//	#          /                                      #
//	#         /                                       #
//	#        /                                        #
//	#       v  x-coord.                               #
//	#                                                 #
//	###################################################

template<class TElement>
SphereGenerator<TElement>::SphereGenerator(int argc, char** argv, int rank, int size)
	: _settings(argc, argv), _rank(rank), _size(size)
{
	if (_settings.layers * 6 != _size) {
		if (_rank == 0) {
			std::cerr << "The number of clusters(" << _settings.layers * 6;
			std::cerr << ") does not accord the number of MPI processes(" << _size << ").\n";
		}
		exit(EXIT_FAILURE);
	}
}

template<class TElement>
void SphereGenerator<TElement>::points(mesh::Coordinates &coordinates)
{
	eslocal cNodes[3];
	SphereUtils<TElement>::clusterNodesCount(_settings, cNodes);

	coordinates.clear();
	coordinates.reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	double step[3];
	for (size_t i = 0; i < 3; i++) {
		step[i] = 1.0 / (cNodes[i] - 1);
	}

	double layerSize = (_settings.outerRadius - _settings.innerRadius) / _settings.layers;
	mesh::Point corner[3];
	corner[0] = mesh::Point(          1,          -1,           1);
	corner[1] = mesh::Point(         -1, corner[0].y, corner[0].z); // x-axis
	corner[2] = mesh::Point(corner[0].x,           1, corner[0].z); // y-axis

	mesh::Point point;
	eslocal index = 0;
	for (eslocal z = 0; z < cNodes[2]; z++) {
		for (eslocal y = 0; y < cNodes[1]; y++) {
			for (eslocal x = 0; x < cNodes[0]; x++) {
				point.z = 1;
				point.y = (corner[0] * (1 - y * step[1]) + corner[2] * y * step[1]).y;
				point.x = (corner[0] * (1 - x * step[0]) + corner[1] * x * step[0]).x;

				point.normalize();
				point *= _settings.innerRadius + layerSize * (_rank / 6) + z * step[2];

				switch (_rank % 6) {
				case 0: { // top
					break;
				}
				case 1: { // bottom
					point.flip();
					break;
				}
				case 2: { // left
					double tmp = point.y;
					point.y = -point.z;
					point.z = tmp;
					break;
				}
				case 3: { // right
					double tmp = point.y;
					point.y = point.z;
					point.z = -tmp;
					break;
				}
				case 4: { // front
					double tmp = point.z;
					point.z = -point.x;
					point.x = tmp;
					break;
				}
				case 5: { // rear
					double tmp = point.z;
					point.z = point.x;
					point.x = -tmp;
					break;
				}
				}

				coordinates.add(point, index, index);
				index++;
			}
		}
	}
}

template<class TElement>
void SphereGenerator<TElement>::elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts)
{
	eslocal cNodes[3];

	SphereUtils<TElement>::clusterNodesCount(_settings, cNodes);

	std::vector<eslocal> indices((2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]) * (2 + TElement::subnodes[2]));

	elements.clear();
	elements.reserve(SphereUtils<TElement>::clusterElementsCount(_settings));
	parts.clear();
	parts.reserve(_settings.subdomainsInCluster[0] * _settings.subdomainsInCluster[1] * _settings.subdomainsInCluster[2] + 1);

	eslocal subdomain[3];
	eslocal element[3];

	eslocal subdomainOffset[3];
	eslocal elementOffset[3];

	parts.push_back(elements.size());
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
												(elementOffset[2] + z) * cNodes[0] * cNodes[1] +
												(elementOffset[1] + y) * cNodes[0] +
												(elementOffset[0] + x);
									}
								}
							}
							e.addElements(elements, &indices[0]);
						}
					}
				}
				parts.push_back(elements.size());
			}
		}
	}
}

template<class TElement>
void SphereGenerator<TElement>::fixPoints(std::vector<eslocal> &fixPoints)
{
	fixPoints.reserve(8 * _settings.subdomainsInCluster[0] * _settings.subdomainsInCluster[1] * _settings.subdomainsInCluster[2]);

	eslocal nodes[3];
	eslocal cNodes[3];
	for (int i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
	}
	SphereUtils<TElement>::clusterNodesCount(_settings, cNodes);

	eslocal offset[3];
	for (eslocal sz = 0; sz < _settings.subdomainsInCluster[2]; sz++) {
		for (eslocal sy = 0; sy < _settings.subdomainsInCluster[1]; sy++) {
			for (eslocal sx = 0; sx < _settings.subdomainsInCluster[0]; sx++) {
				for (int i = 0; i < 8; i++) {
					offset[0] = (i & 1) ? 1 : 0;
					offset[1] = (i & 2) ? 1 : 0;
					offset[2] = (i & 4) ? 1 : 0;
					fixPoints.push_back(
							(sz + offset[2]) * nodes[2] * cNodes[0] * cNodes[1] +
							(sy + offset[1]) * nodes[1] * cNodes[0] +
							(sx + offset[0]) * nodes[0]);
				}
			}
		}
	}
}

template<class TElement>
void SphereGenerator<TElement>::boundaryConditions(mesh::Coordinates &coordinates)
{
	mesh::CoordinatesProperty &dirichlet_x = coordinates.property(mesh::DIRICHLET_X);
	mesh::CoordinatesProperty &dirichlet_y = coordinates.property(mesh::DIRICHLET_Y);
	mesh::CoordinatesProperty &dirichlet_z = coordinates.property(mesh::DIRICHLET_Z);
	mesh::CoordinatesProperty &forces_x = coordinates.property(mesh::FORCES_X);
	mesh::CoordinatesProperty &forces_y = coordinates.property(mesh::FORCES_Y);
	mesh::CoordinatesProperty &forces_z = coordinates.property(mesh::FORCES_Z);

//	eslocal nodes[3];
//	SphereUtils<TElement>::clusterNodesCount(_settings, nodes);
//
//	if (_rank / 6 == 0) {
//		eslocal index = 0;
//		for (eslocal z = 0; z < nodes[2]; z++) {
//			for (eslocal y = 0; y < nodes[1]; y++) {
//				if (_settings.fillCondition[CubeSettings::REAR][mesh::DIRICHLET_X]) {
//					dirichlet_x[index] = _settings.boundaryCondition[CubeSettings::REAR][mesh::DIRICHLET_X];
//				}
//				if (_settings.fillCondition[CubeSettings::REAR][mesh::FORCES_X]) {
//					forces_x[index] = _settings.boundaryCondition[CubeSettings::REAR][mesh::FORCES_X];
//				}
//				if (_settings.fillCondition[CubeSettings::REAR][mesh::DIRICHLET_Y]) {
//					dirichlet_y[index] = _settings.boundaryCondition[CubeSettings::REAR][mesh::DIRICHLET_Y];
//				}
//				if (_settings.fillCondition[CubeSettings::REAR][mesh::FORCES_Y]) {
//					forces_y[index] = _settings.boundaryCondition[CubeSettings::REAR][mesh::FORCES_Y];
//				}
//				if (_settings.fillCondition[CubeSettings::REAR][mesh::DIRICHLET_Z]) {
//					dirichlet_z[index] = _settings.boundaryCondition[CubeSettings::REAR][mesh::DIRICHLET_Z];
//				}
//				if (_settings.fillCondition[CubeSettings::REAR][mesh::FORCES_Z]) {
//					forces_z[index] = _settings.boundaryCondition[CubeSettings::REAR][mesh::FORCES_Z];
//				}
//				index += nodes[0];
//			}
//		}
//	}
}

template <class TElement>
void SphereGenerator<TElement>::corners(mesh::Boundaries &boundaries)
{
	eslocal nodes[3];
	eslocal cNodes[3];
	for (int i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
	}
	SphereUtils<TElement>::clusterNodesCount(_settings, cNodes);

	eslocal step[3];
	for (int i = 0; i < 3; i++) {
		step[i] = _settings.elementsInSubdomain[i] / (_settings.cornerCount + 1);
		step[i] *= TElement::subnodes[i] + 1;
	}
	std::vector<std::vector<size_t> > offsets(3);
	std::vector<size_t> mul(3);

	for (int i = 0; i < 3; i++) {
		for (size_t j = 0; j < _settings.subdomainsInCluster[i]; j++) {
			for (size_t k = 0; k <= _settings.cornerCount / 2; k++) {
				offsets[i].push_back(j * nodes[i] + k * step[i]);
				offsets[i].push_back(j * nodes[i] + nodes[i] - k * step[i]);
			}
			if (_settings.cornerCount % 2 == 1) {
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
					if (!_settings.corners
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] == 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] == 0)
					{
						continue;
					}
					if (!_settings.edges
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] == 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] != 0)
					{
						continue;
					}
					if (!_settings.edges
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] != 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] == 0)
					{
						continue;
					}
					if (!_settings.faces
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] != 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] != 0)
					{
						continue;
					}
					index = i * nodes[d] * mul[d];
					index += offsets[(d + 1) % 3][j] * mul[(d + 1) % 3];
					index += offsets[(d + 2) % 3][k] * mul[(d + 2) % 3];
					boundaries.setCorner(index);
				}
			}
		}
	}
}

template <class TElement>
void SphereGenerator<TElement>::clusterBoundaries(mesh::Boundaries &boundaries)
{
//	esglobal gNodes[3];
//	SphereUtils<TElement>::globalNodesCount(_settings, gNodes);
//	eslocal cNodes[3];
//	SphereUtils<TElement>::clusterNodesCount(_settings, cNodes);
//	boundaries.resize(cNodes[0] * cNodes[1] * cNodes[2]);
//
//	bool border[3];
//	eslocal cIndex = _cluster[0] + _cluster[1] * _settings.clusters[0] + _cluster[2] * _settings.clusters[0] * _settings.clusters[1];
//	esglobal index = 0;
//
//	esglobal cs[3], ce[3];
//	for (eslocal i = 0; i < 3; i++) {
//		cs[i] = (cNodes[i] - 1) * _cluster[i];
//		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
//	}
//
//	for (esglobal z = cs[2]; z <= ce[2]; z++) {
//		border[2] = (z == 0 || z == gNodes[2] - 1) ? false : z % ( cNodes[2] - 1) == 0;
//		for (esglobal y = cs[1]; y <= ce[1]; y++) {
//			border[1] = (y == 0 || y == gNodes[1] - 1) ? false : y % ( cNodes[1] - 1) == 0;
//			for (esglobal x = cs[0]; x <= ce[0]; x++) {
//				border[0] = (x == 0 || x == gNodes[0] - 1) ? false : x % ( cNodes[0] - 1) == 0;
//				for (int i = 0; i < 8; i++) {
//					eslocal tmp = cIndex;
//					if (border[0] && (i & 1)) {
//						tmp += (x == cs[0]) ? -1 : 1;
//					}
//					if (border[1] && (i & 2)) {
//						tmp += ((y == cs[1]) ? -1 : 1) * _settings.clusters[0];
//					}
//					if (border[2] && (i & 4)) {
//						tmp += ((z == cs[2]) ? -1 : 1) * _settings.clusters[0] * _settings.clusters[1];
//					}
//					boundaries[index].insert(tmp);
//				}
//				index++;
//			}
//		}
//	}
}

}

