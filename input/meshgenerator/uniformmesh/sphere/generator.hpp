
#include "generator.h"

namespace espreso {
namespace input {


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

// UP     0
// READ   1
// BOTTOM 2
// FRONT  3
// RIGHT  4
// LEFT   5

static void checkSettings(const SphereSettings &settings)
{
	if (settings.layers * 6 != settings.size) {
		ESINFO(ERROR) << "The number of clusters(" << settings.layers * 6
							<< ") does not accord the number of MPI processes(" << settings.size << ").";
	}
}


template<class TElement>
SphereGenerator<TElement>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<TElement>(mesh, settings), _settings(settings)
{
	checkSettings(_settings);
}

template<class TElement>
void SphereGenerator<TElement>::elementsMaterials(std::vector<Element*> &elements)
{
	// TODO: set materials
}

template<class TElement>
void SphereGenerator<TElement>::points(Coordinates &coordinates, size_t &DOFs)
{
	DOFs = this->_DOFs;
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	coordinates.clear();
	coordinates.reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	double step[3];
	for (size_t i = 0; i < 3; i++) {
		step[i] = 1.0 / (cNodes[i] - 1);
	}

	double layerSize = (_settings.outerRadius - _settings.innerRadius) / _settings.layers;
	Point corner[3];
	corner[0] = Point(          1,          -1,           1);
	corner[1] = Point(         -1, corner[0].y, corner[0].z); // x-axis
	corner[2] = Point(corner[0].x,           1, corner[0].z); // y-axis

	eslocal surface = SphereUtils<TElement>::surfaceNodesCount(_settings);
	eslocal ring = SphereUtils<TElement>::ringNodesCount(_settings);

	Point point;
	eslocal index = 0;
	esglobal gIndex;
	for (eslocal z = 0; z < cNodes[2]; z++) {
		for (eslocal y = 0; y < cNodes[1]; y++) {
			for (eslocal x = cNodes[0] - 1; x >= 0; x--) {
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

				gIndex += _settings.index / 6 * surface * (cNodes[2] - 1);
				coordinates.add(point, index++, gIndex);
			}
		}
	}
}


template<class TElement>
void SphereGenerator<TElement>::boundaryConditions(Coordinates &coordinates)
{
	if (_settings.index > 5) {
		return;
	}
	CoordinatesProperty &dirichlet_x = coordinates.property(DIRICHLET_X);
	CoordinatesProperty &dirichlet_y = coordinates.property(DIRICHLET_Y);
	CoordinatesProperty &dirichlet_z = coordinates.property(DIRICHLET_Z);

	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	for (eslocal i = 0; i < cNodes[0] * cNodes[1]; i++) {
		dirichlet_x[i] = 0;
		dirichlet_y[i] = 0;
		dirichlet_z[i] = 0;
	}
}


template <class TElement>
void SphereGenerator<TElement>::clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours)
{
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	boundaries.resize(cNodes[0] * cNodes[1] * cNodes[2]);

	std::vector<eslocal> clusterMap(15);
	for (int i = -1; i < 2; i++) {
		switch (_settings.index % 6) {
		case 0: case 1: case 2: case 3:
			clusterMap[5 + i * 5] = i * 6 + _settings.index / 6 * 6 + (_settings.index % 6 % 4 + 3) % 4;
			clusterMap[6 + i * 5] = i * 6 + _settings.index / 6 * 6 + 5;
			clusterMap[7 + i * 5] = i * 6 + _settings.index;
			clusterMap[8 + i * 5] = i * 6 + _settings.index / 6 * 6 + 4;
			clusterMap[9 + i * 5] = i * 6 + _settings.index / 6 * 6 + (_settings.index % 6 % 4 + 1) % 4;
			break;
		case 4:
			clusterMap[5 + i * 5] = i * 6 + _settings.index / 6 * 6 + 3;
			clusterMap[6 + i * 5] = i * 6 + _settings.index / 6 * 6 + 0;
			clusterMap[7 + i * 5] = i * 6 + _settings.index;
			clusterMap[8 + i * 5] = i * 6 + _settings.index / 6 * 6 + 2;
			clusterMap[9 + i * 5] = i * 6 + _settings.index / 6 * 6 + 1;
			break;
		case 5:
			clusterMap[5 + i * 5] = i * 6 + _settings.index / 6 * 6 + 3;
			clusterMap[6 + i * 5] = i * 6 + _settings.index / 6 * 6 + 2;
			clusterMap[7 + i * 5] = i * 6 + _settings.index;
			clusterMap[8 + i * 5] = i * 6 + _settings.index / 6 * 6 + 0;
			clusterMap[9 + i * 5] = i * 6 + _settings.index / 6 * 6 + 1;
			break;
		}
	}

	std::vector<eslocal> clusters(clusterMap);
	std::sort(clusters.begin(), clusters.end());

	for (size_t i = 0; i < 15; i++) {
		if (0 <= clusters[i] && clusters[i] < _settings.size) {
			size_t index = std::find(clusterMap.begin(), clusterMap.end(), clusters[i]) - clusterMap.begin();
			eslocal sx = index % 5 == 3 ? cNodes[0] - 1 : 0,
					ex = index % 5 != 1 ? cNodes[0] - 1 : 0,
					sy = index % 5 == 4 ? cNodes[1] - 1 : 0,
					ey = index % 5 != 0 ? cNodes[1] - 1 : 0,
					sz = index >= 10    ? cNodes[2] - 1 : 0,
					ez = index >= 5     ? cNodes[2] - 1 : 0;

			for (eslocal z = sz; z <= ez; z++) {
				for (eslocal y = sy; y <= ey; y++) {
					for (eslocal x = sx; x <= ex; x++) {
						boundaries[z * cNodes[0] * cNodes[1] + y * cNodes[0] + x].push_back(clusters[i]);
					}
				}
			}
			if (clusters[i] != config::MPIrank) {
				neighbours.push_back(clusters[i]);
			}
		}
	}
}

}
}

