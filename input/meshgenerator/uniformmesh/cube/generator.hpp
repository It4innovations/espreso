
#include "generator.h"

namespace esinput {


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

static void setCluster(size_t cluster[], const CubeSettings &settings)
{
	if (settings.clusters[0] * settings.clusters[1] * settings.clusters[2] != settings.size) {
		if (settings.index == 0) {
			std::cerr << "The number of clusters(" << settings.clusters[0] * settings.clusters[1] * settings.clusters[2];
			std::cerr << ") does not accord the number of MPI processes(" << settings.size << ").\n";
		}
		exit(EXIT_FAILURE);
	}
	eslocal index = 0, i = 0;
	for (size_t z = 0; z < settings.clusters[2]; z++) {
		for (size_t y = 0; y < settings.clusters[1]; y++) {
			for (size_t x = 0; x < settings.clusters[0]; x++) {
				if (settings.index == index++) {
					cluster[0] = x;
					cluster[1] = y;
					cluster[2] = z;
					return;
				}
			}
		}
	}
}

template<class TElement>
CubeGenerator<TElement>::CubeGenerator(int argc, char** argv, size_t index, size_t size)
	: UniformGenerator<TElement>(argc, argv, index, size), _settings(argc, argv, index, size)
{
	setCluster(_cluster, _settings);
}

template<class TElement>
CubeGenerator<TElement>::CubeGenerator(const CubeSettings &settings)
	: UniformGenerator<TElement>(settings), _settings(settings)
{
	setCluster(_cluster, _settings);
}

template<class TElement>
void CubeGenerator<TElement>::elementsMaterials(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts)
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
									elements[counter++]->setParam(mesh::Element::MATERIAL, material);
								}
							}
						}
					}

				}
			}
	}

}


template<class TElement>
void CubeGenerator<TElement>::points(mesh::Coordinates &coordinates)
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
					mesh::Point(x * step[0], y * step[1], z * step[2]),
					(z - cs[2]) * cNodes[0] * cNodes[1] + (y - cs[1]) * cNodes[0] + (x - cs[0]),
					z * gNodes[0] * gNodes[1] + y * gNodes[0] + x
				);
			}
		}
	}
}

template<class TElement>
void CubeGenerator<TElement>::boundaryConditions(mesh::Coordinates &coordinates)
{
	mesh::CoordinatesProperty &dirichlet_x = coordinates.property(mesh::DIRICHLET_X);
	mesh::CoordinatesProperty &dirichlet_y = coordinates.property(mesh::DIRICHLET_Y);
	mesh::CoordinatesProperty &dirichlet_z = coordinates.property(mesh::DIRICHLET_Z);
	mesh::CoordinatesProperty &forces_x = coordinates.property(mesh::FORCES_X);
	mesh::CoordinatesProperty &forces_y = coordinates.property(mesh::FORCES_Y);
	mesh::CoordinatesProperty &forces_z = coordinates.property(mesh::FORCES_Z);

	eslocal nodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, nodes);

	if (_cluster[0] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal y = 0; y < nodes[1]; y++) {
				if (_settings.fillCondition[CubeSettings::REAR].find(mesh::DIRICHLET_X)->second) {
					dirichlet_x[index] = _settings.boundaryCondition[CubeSettings::REAR].find(mesh::DIRICHLET_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::REAR].find(mesh::FORCES_X)->second) {
					forces_x[index] = _settings.boundaryCondition[CubeSettings::REAR].find(mesh::FORCES_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::REAR].find(mesh::DIRICHLET_Y)->second) {
					dirichlet_y[index] = _settings.boundaryCondition[CubeSettings::REAR].find(mesh::DIRICHLET_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::REAR].find(mesh::FORCES_Y)->second) {
					forces_y[index] = _settings.boundaryCondition[CubeSettings::REAR].find(mesh::FORCES_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::REAR].find(mesh::DIRICHLET_Z)->second) {
					dirichlet_z[index] = _settings.boundaryCondition[CubeSettings::REAR].find(mesh::DIRICHLET_Z)->second;
				}
				if (_settings.fillCondition[CubeSettings::REAR].find(mesh::FORCES_Z)->second) {
					forces_z[index] = _settings.boundaryCondition[CubeSettings::REAR].find(mesh::FORCES_Z)->second;
				}
				index += nodes[0];
			}
		}
	}

	if (_cluster[0] == _settings.clusters[0] - 1) {
		eslocal index = nodes[0] - 1;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal y = 0; y < nodes[1]; y++) {
				if (_settings.fillCondition[CubeSettings::FRONT].find(mesh::DIRICHLET_X)->second) {
					dirichlet_x[index] = _settings.boundaryCondition[CubeSettings::FRONT].find(mesh::DIRICHLET_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::FRONT].find(mesh::FORCES_X)->second) {
					forces_x[index] = _settings.boundaryCondition[CubeSettings::FRONT].find(mesh::FORCES_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::FRONT].find(mesh::DIRICHLET_Y)->second) {
					dirichlet_y[index] = _settings.boundaryCondition[CubeSettings::FRONT].find(mesh::DIRICHLET_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::FRONT].find(mesh::FORCES_Y)->second) {
					forces_y[index] = _settings.boundaryCondition[CubeSettings::FRONT].find(mesh::FORCES_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::FRONT].find(mesh::DIRICHLET_Z)->second) {
					dirichlet_z[index] = _settings.boundaryCondition[CubeSettings::FRONT].find(mesh::DIRICHLET_Z)->second;
				}
				if (_settings.fillCondition[CubeSettings::FRONT].find(mesh::FORCES_Z)->second) {
					forces_z[index] = _settings.boundaryCondition[CubeSettings::FRONT].find(mesh::FORCES_Z)->second;
				}
				index += nodes[0];
			}
		}
	}

	if (_cluster[1] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (_settings.fillCondition[CubeSettings::LEFT].find(mesh::DIRICHLET_X)->second) {
					dirichlet_x[index] = _settings.boundaryCondition[CubeSettings::LEFT].find(mesh::DIRICHLET_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::LEFT].find(mesh::FORCES_X)->second) {
					forces_x[index] = _settings.boundaryCondition[CubeSettings::LEFT].find(mesh::FORCES_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::LEFT].find(mesh::DIRICHLET_Y)->second) {
					dirichlet_y[index] = _settings.boundaryCondition[CubeSettings::LEFT].find(mesh::DIRICHLET_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::LEFT].find(mesh::FORCES_Y)->second) {
					forces_y[index] = _settings.boundaryCondition[CubeSettings::LEFT].find(mesh::FORCES_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::LEFT].find(mesh::DIRICHLET_Z)->second) {
					dirichlet_z[index] = _settings.boundaryCondition[CubeSettings::LEFT].find(mesh::DIRICHLET_Z)->second;
				}
				if (_settings.fillCondition[CubeSettings::LEFT].find(mesh::FORCES_Z)->second) {
					forces_z[index] = _settings.boundaryCondition[CubeSettings::LEFT].find(mesh::FORCES_Z)->second;
				}
				index++;
			}
			index = (z + 1) * nodes[1] * nodes[0];
		}
	}

	if (_cluster[1] == _settings.clusters[1] - 1) {
		eslocal index = nodes[1] * nodes[0] - 1;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (_settings.fillCondition[CubeSettings::RIGHT].find(mesh::DIRICHLET_X)->second) {
					dirichlet_x[index] = _settings.boundaryCondition[CubeSettings::RIGHT].find(mesh::DIRICHLET_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::RIGHT].find(mesh::FORCES_X)->second) {
					forces_x[index] = _settings.boundaryCondition[CubeSettings::RIGHT].find(mesh::FORCES_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::RIGHT].find(mesh::DIRICHLET_Y)->second) {
					dirichlet_y[index] = _settings.boundaryCondition[CubeSettings::RIGHT].find(mesh::DIRICHLET_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::RIGHT].find(mesh::FORCES_Y)->second) {
					forces_y[index] = _settings.boundaryCondition[CubeSettings::RIGHT].find(mesh::FORCES_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::RIGHT].find(mesh::DIRICHLET_Z)->second) {
					dirichlet_z[index] = _settings.boundaryCondition[CubeSettings::RIGHT].find(mesh::DIRICHLET_Z)->second;
				}
				if (_settings.fillCondition[CubeSettings::RIGHT].find(mesh::FORCES_Z)->second) {
					forces_z[index] = _settings.boundaryCondition[CubeSettings::RIGHT].find(mesh::FORCES_Z)->second;
				}
				index++;
			}
			index = (z + 2) * nodes[1] * nodes[0] - 1;
		}
	}

	if (_cluster[2] == 0) {
		eslocal index = 0;
		for (eslocal y = 0; y < nodes[1]; y++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (_settings.fillCondition[CubeSettings::BOTTOM].find(mesh::DIRICHLET_X)->second) {
					dirichlet_x[index] = _settings.boundaryCondition[CubeSettings::BOTTOM].find(mesh::DIRICHLET_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::BOTTOM].find(mesh::FORCES_X)->second) {
					forces_x[index] = _settings.boundaryCondition[CubeSettings::BOTTOM].find(mesh::FORCES_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::BOTTOM].find(mesh::DIRICHLET_Y)->second) {
					dirichlet_y[index] = _settings.boundaryCondition[CubeSettings::BOTTOM].find(mesh::DIRICHLET_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::BOTTOM].find(mesh::FORCES_Y)->second) {
					forces_y[index] = _settings.boundaryCondition[CubeSettings::BOTTOM].find(mesh::FORCES_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::BOTTOM].find(mesh::DIRICHLET_Z)->second) {
					dirichlet_z[index] = _settings.boundaryCondition[CubeSettings::BOTTOM].find(mesh::DIRICHLET_Z)->second;
				}
				if (_settings.fillCondition[CubeSettings::BOTTOM].find(mesh::FORCES_Z)->second) {
					forces_z[index] = _settings.boundaryCondition[CubeSettings::BOTTOM].find(mesh::FORCES_Z)->second;
				}
				index++;
			}
		}
	}

	if (_cluster[2] == _settings.clusters[2] - 1) {
		eslocal index = nodes[0] * nodes[1] * (nodes[2] - 1);
		for (eslocal y = 0; y < nodes[1]; y++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (_settings.fillCondition[CubeSettings::TOP].find(mesh::DIRICHLET_X)->second) {
					dirichlet_x[index] = _settings.boundaryCondition[CubeSettings::TOP].find(mesh::DIRICHLET_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::TOP].find(mesh::FORCES_X)->second) {
					forces_x[index] = _settings.boundaryCondition[CubeSettings::TOP].find(mesh::FORCES_X)->second;
				}
				if (_settings.fillCondition[CubeSettings::TOP].find(mesh::DIRICHLET_Y)->second) {
					dirichlet_y[index] = _settings.boundaryCondition[CubeSettings::TOP].find(mesh::DIRICHLET_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::TOP].find(mesh::FORCES_Y)->second) {
					forces_y[index] = _settings.boundaryCondition[CubeSettings::TOP].find(mesh::FORCES_Y)->second;
				}
				if (_settings.fillCondition[CubeSettings::TOP].find(mesh::DIRICHLET_Z)->second) {
					dirichlet_z[index] = _settings.boundaryCondition[CubeSettings::TOP].find(mesh::DIRICHLET_Z)->second;
				}
				if (_settings.fillCondition[CubeSettings::TOP].find(mesh::FORCES_Z)->second) {
					forces_z[index] = _settings.boundaryCondition[CubeSettings::TOP].find(mesh::FORCES_Z)->second;
				}
				index++;
			}
		}
	}
}


template <class TElement>
void CubeGenerator<TElement>::clusterBoundaries(mesh::Boundaries &boundaries)
{
	esglobal gNodes[3];
	CubeUtils<TElement>::globalNodesCount(_settings, gNodes);
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	boundaries.resize(cNodes[0] * cNodes[1] * cNodes[2]);

	bool border[3];
	eslocal cIndex = _cluster[0] + _cluster[1] * _settings.clusters[0] + _cluster[2] * _settings.clusters[0] * _settings.clusters[1];
	esglobal index = 0;

	esglobal cs[3], ce[3];
	for (eslocal i = 0; i < 3; i++) {
		cs[i] = (cNodes[i] - 1) * _cluster[i];
		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
	}

	for (esglobal z = cs[2]; z <= ce[2]; z++) {
		border[2] = (z == 0 || z == gNodes[2] - 1) ? false : z % ( cNodes[2] - 1) == 0;
		for (esglobal y = cs[1]; y <= ce[1]; y++) {
			border[1] = (y == 0 || y == gNodes[1] - 1) ? false : y % ( cNodes[1] - 1) == 0;
			for (esglobal x = cs[0]; x <= ce[0]; x++) {
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

}

