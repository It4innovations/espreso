
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

static void checkSettings(const SphereSettings &settings)
{
	if (settings.layers * 6 != settings.size) {
		if (settings.index == 0) {
			std::cerr << "The number of clusters(" << settings.layers * 6;
			std::cerr << ") does not accord the number of MPI processes(" << settings.size << ").\n";
		}
		exit(EXIT_FAILURE);
	}
}

template<class TElement>
SphereGenerator<TElement>::SphereGenerator(int argc, char** argv, size_t index, size_t size)
	: UniformGenerator<TElement>(argc, argv, index, size), _settings(argc, argv, index, size)
{
	checkSettings(_settings);
}

template<class TElement>
SphereGenerator<TElement>::SphereGenerator(const SphereSettings &settings)
	: UniformGenerator<TElement>(_settings), _settings(settings)
{
	checkSettings(_settings);
}

template<class TElement>
void SphereGenerator<TElement>::elementsMaterials(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts)
{
	// TODO: set materials
}

template<class TElement>
void SphereGenerator<TElement>::points(mesh::Coordinates &coordinates)
{
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

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

	eslocal surface = SphereUtils<TElement>::surfaceNodesCount(_settings);
	eslocal ring = SphereUtils<TElement>::ringNodesCount(_settings);

	mesh::Point point;
	eslocal index = 0;
	esglobal gIndex;
	for (eslocal z = 0; z < cNodes[2]; z++) {
		for (eslocal x = 0; x < cNodes[0]; x++) {
			for (eslocal y = 0; y < cNodes[1]; y++) {
				point.z = 1;
				point.y = (corner[0] * (1 - y * step[1]) + corner[2] * y * step[1]).y;
				point.x = (corner[0] * (1 - x * step[0]) + corner[1] * x * step[0]).x;

				point.normalize();
				point *= _settings.innerRadius + layerSize * ((_settings.index / 6) + z * step[2]);

				switch (_settings.index % 6) {
				case 0: { // top
					gIndex = z * surface + x * ring + y;
					break;
				}
				case 1: { // right
					gIndex = z * surface + x * ring + y + cNodes[1] - 1;
					double tmp = point.y;
					point.y = point.z;
					point.z = -tmp;
					break;
				}
				case 2: { // bottom
					gIndex = z * surface + x * ring + y + 2 * cNodes[1] - 2;
					point.z = -point.z;
					point.y = -point.y;
					break;
				}
				case 3: { // left
					gIndex = z * surface + x * ring + (y + 3 * cNodes[1] - 3) % ring;
					double tmp = point.y;
					point.y = -point.z;
					point.z = tmp;
					break;
				}
				case 4: { // front
					gIndex = z * surface + cNodes[0] * ring + (x - 1) * (cNodes[1] - 2) + y - 1;
					if (y == 0) {
						gIndex = z * surface + (x + 3 * cNodes[1] - 3) % ring;
					}
					if (y == cNodes[1] - 1) {
						gIndex = z * surface + cNodes[0] - 1 - x + cNodes[1] - 1;
					}
					if (x == 0) {
						gIndex = z * surface + x * ring + cNodes[1] - 1 - y + 2 * cNodes[1] - 2;
					}
					if (x == cNodes[0] - 1) {
						gIndex = z * surface + y;
					}
					double tmp = point.z;
					point.z = -point.x;
					point.x = tmp;
					break;
				}
				case 5: { // rear
					gIndex = z * surface + cNodes[0] * ring + (x - 1) * (cNodes[1] - 2) + y - 1 + (cNodes[1] - 2) * (cNodes[0] - 2);
					if (y == 0) {
						gIndex = z * surface + (cNodes[0] - 1) * ring + (cNodes[1] - 1 - x + 3 * cNodes[1] - 3) % ring;
					}
					if (y == cNodes[1] - 1) {
						gIndex = z * surface + (cNodes[0] - 1) * ring + x + cNodes[1] - 1;
					}
					if (x == 0) {
						gIndex = z * surface + (cNodes[0] - 1) * ring + y;
					}
					if (x == cNodes[0] - 1) {
						gIndex = z * surface + x * ring + cNodes[1] - 1 - y + 2 * cNodes[1] - 2;
					}
					double tmp = point.z;
					point.z = point.x;
					point.x = -tmp;
					break;
				}
				}

				coordinates.add(point, index++, gIndex);
			}
		}
	}
}


template<class TElement>
void SphereGenerator<TElement>::boundaryConditions(mesh::Coordinates &coordinates)
{
	if (_settings.index > 5) {
		return;
	}
	mesh::CoordinatesProperty &dirichlet_x = coordinates.property(mesh::DIRICHLET_X);
	mesh::CoordinatesProperty &dirichlet_y = coordinates.property(mesh::DIRICHLET_Y);
	mesh::CoordinatesProperty &dirichlet_z = coordinates.property(mesh::DIRICHLET_Z);

	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	for (eslocal i = 0; i < cNodes[0] * cNodes[1]; i++) {
		dirichlet_x[i] = 0;
		dirichlet_y[i] = 0;
		dirichlet_z[i] = 0;
	}
}


template <class TElement>
void SphereGenerator<TElement>::clusterBoundaries(mesh::Boundaries &boundaries)
{
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	boundaries.resize(cNodes[0] * cNodes[1] * cNodes[2]);

	for (size_t i = 0; i < boundaries.size(); i++) {
		boundaries[i].insert(_settings.index);
	}

	if (_settings.index / 6 > 0) {
		for (eslocal i = 0; i < cNodes[1] * cNodes[0]; i++) {
			boundaries[i].insert(_settings.index - 6);
		}
	}

	if (_settings.index / 6 < _settings.layers - 1) {
		for (eslocal i = 0; i < cNodes[1] * cNodes[0]; i++) {
			boundaries[i + (cNodes[2] - 1) * cNodes[1] * cNodes[0]].insert(_settings.index + 6);
		}
	}

	eslocal index;
	if  (_settings.index % 6 < 4) {
		for (eslocal z = 0; z < cNodes[2]; z++) {
			index = z * cNodes[0] * cNodes[1];
			for (eslocal y = 0; y < cNodes[1]; y++) {
				boundaries[index + y].insert(6 * (_settings.index / 6) + 4);
				boundaries[index + y + (cNodes[0] - 1) * cNodes[1]].insert(6 * (_settings.index / 6) + 5);
			}
			for (eslocal x = 0; x < cNodes[0]; x++) {
				boundaries[index + x * cNodes[1]].insert(6 * (_settings.index / 6) + (_settings.index + 3) % 4);
				boundaries[index + (x + 1) * cNodes[1] - 1].insert(6 * (_settings.index / 6) + (_settings.index + 1) % 4);
			}
		}
	} else if (_settings.index % 6 == 4) {
		for (eslocal z = 0; z < cNodes[2]; z++) {
			index = z * cNodes[0] * cNodes[1];
			for (eslocal y = 0; y < cNodes[1]; y++) {
				boundaries[index + y].insert(6 * (_settings.index / 6) + 2);
				boundaries[index + y + (cNodes[0] - 1) * cNodes[1]].insert(6 * (_settings.index / 6));
			}
			for (eslocal x = 0; x < cNodes[0]; x++) {
				boundaries[index + x * cNodes[1]].insert(6 * (_settings.index / 6) + 3);
				boundaries[index + (x + 1) * cNodes[1] - 1].insert(6 * (_settings.index / 6) + 1);
			}
		}
	} else {
		for (eslocal z = 0; z < cNodes[2]; z++) {
			index = z * cNodes[0] * cNodes[1];
			for (eslocal y = 0; y < cNodes[1]; y++) {
				boundaries[index + y].insert(6 * (_settings.index / 6));
				boundaries[index + y + (cNodes[0] - 1) * cNodes[1]].insert(6 * (_settings.index / 6) + 2);
			}
			for (eslocal x = 0; x < cNodes[0]; x++) {
				boundaries[index + x * cNodes[1]].insert(6 * (_settings.index / 6) + 3);
				boundaries[index + (x + 1) * cNodes[1] - 1].insert(6 * (_settings.index / 6) + 1);
			}
		}
	}
}

}

