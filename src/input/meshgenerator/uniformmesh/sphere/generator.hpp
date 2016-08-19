
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
// RIGHT  1
// BOTTOM 2
// LEFT   3
// FRONT  4
// REAR   5

template<class TElement>
SphereGenerator<TElement>::SphereGenerator(Mesh &mesh, const SphereSettings &settings)
	: UniformGenerator<TElement>(mesh, settings), _settings(settings)
{
	ESINFO(GLOBAL_ERROR) << "Sphere generator does not support the selected element type.";
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

	// Idea: merge clusters on each side together for computation of global indices
	SphereSettings merged = _settings;
	for (size_t i = 0; i < 2; i++) {
		merged.subdomainsInCluster[i] *= _settings.grid;
	}
	merged.subdomainsInCluster[2] *= _settings.layers;
	eslocal sNodes[3]; // side nodes counters
	UniformUtils<TElement>::clusterNodesCount(merged, sNodes);

	coordinates.clear();
	coordinates.reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	double step[3];
	for (size_t i = 0; i < 2; i++) {
		step[i] = 1.0 / (sNodes[i] - 1);
	}
	step[2] = (_settings.outerRadius - _settings.innerRadius) / (sNodes[2] - 1);

	Point corner[3];
	corner[0] = Point(          1,          -1,           1);
	corner[1] = Point(         -1, corner[0].y, corner[0].z); // x-axis
	corner[2] = Point(corner[0].x,           1, corner[0].z); // y-axis

	eslocal surface = SphereUtils<TElement>::surfaceNodesCount(merged);
	eslocal ring = SphereUtils<TElement>::ringNodesCount(merged);

	Point point;
	eslocal index = 0;
	esglobal gIndex;

	eslocal eNodes[2] = { (cNodes[0] - 1) * _settings.grid + 1, (cNodes[1] - 1) * _settings.grid + 1 };

	for (eslocal z = (cNodes[2] - 1) * _cluster[2]; z < (cNodes[2] - 1) * _cluster[2] + cNodes[2]; z++) {
		for (eslocal y = (cNodes[1] - 1) * _cluster[1]; y < (cNodes[1] - 1) * _cluster[1] + cNodes[1]; y++) {
			for (eslocal x = (cNodes[0] - 1) * _cluster[0]; x < (cNodes[0] - 1) * _cluster[0] + cNodes[0]; x++) {
				point.z = 1;
				point.y = (corner[0] * (1 - y * step[1]) + corner[2] * y * step[1]).y;
				point.x = (corner[0] * x * step[0] + corner[1] * (1 - x * step[0])).x;

				point.normalize();
				point *= _settings.innerRadius + z * step[2];

				switch (_side) {
				case 0: { // top
					gIndex = z * surface + x * ring + y;
					break;
				}
				case 1: { // right
					gIndex = z * surface + x * ring + y + eNodes[1] - 1;
					double tmp = point.y;
					point.y = point.z;
					point.z = -tmp;
					break;
				}
				case 2: { // bottom
					gIndex = z * surface + x * ring + y + 2 * eNodes[1] - 2;
					point.z = -point.z;
					point.y = -point.y;
					break;
				}
				case 3: { // left
					gIndex = z * surface + x * ring + (y + 3 * eNodes[1] - 3) % ring;
					double tmp = point.y;
					point.y = -point.z;
					point.z = tmp;
					break;
				}
				case 4: { // front
					gIndex = z * surface + eNodes[0] * ring + (x - 1) * (eNodes[1] - 2) + y - 1;
					if (y == 0) {
						gIndex = z * surface + (eNodes[0] - 1) * ring + 4 * eNodes[1] - 4 - x;
					}
					if (y == eNodes[1] - 1) {
						gIndex = z * surface + (eNodes[0] - 1) * ring + 1 * eNodes[1] - 1 + x;
					}
					if (x == 0) {
						gIndex = z * surface + (eNodes[0] - 1) * ring                     + y;
					}
					if (x == eNodes[0] - 1) {
						gIndex = z * surface + (eNodes[0] - 1) * ring + 3 * eNodes[1] - 3 - y;
					}
					double tmp = point.z;
					point.z = -point.x;
					point.x = tmp;
					break;
				}
				case 5: { // rear
					gIndex = z * surface + eNodes[0] * ring + (x - 1) * (eNodes[1] - 2) + y - 1 + (eNodes[1] - 2) * (eNodes[0] - 2);
					if (y == 0) {
						gIndex = z * surface + 3 * eNodes[1] - 3 + x;
					}
					if (y == eNodes[1] - 1) {
						gIndex = z * surface + 2 * eNodes[1] - 2 - x;
					}
					if (x == 0) {
						gIndex = z * surface + 3 * eNodes[1] - 3 - y;
					}
					if (x == eNodes[0] - 1) {
						gIndex = z * surface                     + y;
					}
					double tmp = point.z;
					point.z = point.x;
					point.x = -tmp;
					break;
				}
				}

				//gIndex += _settings.index / 6 * surface * (cNodes[2] - 1);
				coordinates.add(point, index++, gIndex);
			}
		}
	}
}


/*
 *
 *
 *               | z
 *               |
 *               |
 *               |
 *               -------------y
 *              /
 *             /
 *            /
 *           /x
 */
template <class TElement>
void SphereGenerator<TElement>::clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours)
{
	eslocal line   = _settings.grid;
	eslocal square = _settings.grid * line;
	eslocal cube   = _settings.layers * square;

	std::vector<int> neighboursMap(27); // 3 x 3 x 3

	for (int i = _cluster[2] ? -1 : 0; i < 2; i++) {
		neighboursMap[13 + 9 * i] = _settings.index + i * square;
		neighboursMap[12 + 9 * i] = neighboursMap[13 + 9 * i] - 1;
		neighboursMap[14 + 9 * i] = neighboursMap[13 + 9 * i] + 1;
		neighboursMap[10 + 9 * i] = neighboursMap[13 + 9 * i] - line;
		neighboursMap[16 + 9 * i] = neighboursMap[13 + 9 * i] + line;
		neighboursMap[ 9 + 9 * i] = neighboursMap[10 + 9 * i] - 1;
		neighboursMap[11 + 9 * i] = neighboursMap[10 + 9 * i] + 1;
		neighboursMap[15 + 9 * i] = neighboursMap[16 + 9 * i] - 1;
		neighboursMap[17 + 9 * i] = neighboursMap[16 + 9 * i] + 1;
		if (_cluster[0] == 0) {
			switch (_side) {
			case 0:
				neighboursMap[ 9 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + line * (_cluster[1] + 0) - 1;
				neighboursMap[12 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + line * (_cluster[1] + 1) - 1;
				neighboursMap[15 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + line * (_cluster[1] + 2) - 1;
				break;
			case 1:
				neighboursMap[ 9 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + square - line + _settings.grid - _cluster[1] - 0;
				neighboursMap[12 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + square - line + _settings.grid - _cluster[1] - 1;
				neighboursMap[15 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + square - line + _settings.grid - _cluster[1] - 2;
				break;
			case 2:
				neighboursMap[ 9 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 0);
				neighboursMap[12 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 1);
				neighboursMap[15 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 2);
				break;
			case 3:
				neighboursMap[ 9 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + _cluster[1] - 1;
				neighboursMap[12 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + _cluster[1] + 0;
				neighboursMap[15 + 9 * i] = 5 * cube + (_cluster[2] + i) * square + _cluster[1] + 1;
				break;
			case 4:
				neighboursMap[ 9 + 9 * i] = (_cluster[2] + i) * square + line * (_cluster[1] + 0) - 1;
				neighboursMap[12 + 9 * i] = (_cluster[2] + i) * square + line * (_cluster[1] + 1) - 1;
				neighboursMap[15 + 9 * i] = (_cluster[2] + i) * square + line * (_cluster[1] + 2) - 1;
				break;
			case 5:
				neighboursMap[ 9 + 9 * i] = 2 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 0);
				neighboursMap[12 + 9 * i] = 2 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 1);
				neighboursMap[15 + 9 * i] = 2 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 2);
				break;
			}
		}
		if (_cluster[0] + 1 == _settings.grid) {
			switch (_side) {
			case 0:
				neighboursMap[11 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + line * (_cluster[1] - 1);
				neighboursMap[14 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + line * (_cluster[1] + 0);
				neighboursMap[17 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + line * (_cluster[1] + 1);
				break;
			case 1:
				neighboursMap[11 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + square - line + _cluster[1] - 1;
				neighboursMap[14 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + square - line + _cluster[1] + 0;
				neighboursMap[17 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + square - line + _cluster[1] + 1;
				break;
			case 2:
				neighboursMap[11 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] + 1) - 1;
				neighboursMap[14 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] + 0) - 1;
				neighboursMap[17 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 1) - 1;
				break;
			case 3:
				neighboursMap[11 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + _settings.grid - _cluster[1] - 0;
				neighboursMap[14 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + _settings.grid - _cluster[1] - 1;
				neighboursMap[17 + 9 * i] = 4 * cube + (_cluster[2] + i) * square + _settings.grid - _cluster[1] - 2;
				break;
			case 4:
				neighboursMap[11 + 9 * i] = 2 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] + 1) - 1;
				neighboursMap[14 + 9 * i] = 2 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] + 0) - 1;
				neighboursMap[17 + 9 * i] = 2 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[1] - 1) - 1;
				break;
			case 5:
				neighboursMap[11 + 9 * i] = (_cluster[2] + i) * square + line * (_cluster[1] - 1);
				neighboursMap[14 + 9 * i] = (_cluster[2] + i) * square + line * (_cluster[1] + 0);
				neighboursMap[17 + 9 * i] = (_cluster[2] + i) * square + line * (_cluster[1] + 1);
				break;
			}
		}
		if (_cluster[1] == 0) {
			switch (_side) {
			case 0: case 1: case 2: case 3:
				neighboursMap[10 + 9 * i] = ((_side + 3) % 4) * cube + (_cluster[2] + i) * square + square - line + _cluster[0];
				neighboursMap[ 9 + 9 * i] = neighboursMap[10 + 9 * i] - 1;
				neighboursMap[11 + 9 * i] = neighboursMap[10 + 9 * i] + 1;
				break;
			case 4:
				neighboursMap[10 + 9 * i] = 3 * cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[0]) - 1;
				neighboursMap[ 9 + 9 * i] = neighboursMap[10 + 9 * i] + line;
				neighboursMap[11 + 9 * i] = neighboursMap[10 + 9 * i] - line;
				break;
			case 5:
				neighboursMap[10 + 9 * i] = 3 * cube + (_cluster[2] + i) * square + line *  _cluster[0];
				neighboursMap[ 9 + 9 * i] = neighboursMap[10 + 9 * i] - line;
				neighboursMap[11 + 9 * i] = neighboursMap[10 + 9 * i] + line;
				break;
			}
		}
		if (_cluster[1] + 1 == _settings.grid) {
			switch (_side) {
			case 0: case 1: case 2: case 3:
				neighboursMap[16 + 9 * i] = ((_side + 1) % 4) * cube + (_cluster[2] + i) * square + _cluster[0];
				neighboursMap[15 + 9 * i] = neighboursMap[16 + 9 * i] - 1;
				neighboursMap[17 + 9 * i] = neighboursMap[16 + 9 * i] + 1;
				break;
			case 4:
				neighboursMap[16 + 9 * i] = cube + (_cluster[2] + i) * square + line * (_cluster[0] + 1) - 1;
				neighboursMap[15 + 9 * i] = neighboursMap[16 + 9 * i] - line;
				neighboursMap[17 + 9 * i] = neighboursMap[16 + 9 * i] + line;
				break;
			case 5:
				neighboursMap[16 + 9 * i] = cube + (_cluster[2] + i) * square + line * (_settings.grid - _cluster[0] - 1);
				neighboursMap[15 + 9 * i] = neighboursMap[16 + 9 * i] + line;
				neighboursMap[17 + 9 * i] = neighboursMap[16 + 9 * i] - line;
				break;
			}
		}
		if (_cluster[0] == 0 && _cluster[1] == 0) {
			neighboursMap[9 + 9 * i] = -1;
		}
		if (_cluster[0] == 0 && _cluster[1] + 1 == _settings.grid) {
			neighboursMap[15 + 9 * i] = -1;
		}
		if (_cluster[0] + 1 == _settings.grid && _cluster[1] == 0) {
			neighboursMap[11 + 9 * i] = -1;
		}
		if (_cluster[0] + 1 == _settings.grid && _cluster[1] + 1 == _settings.grid) {
			neighboursMap[17 + 9 * i] = -1;
		}
	}


	if (_cluster[2] == 0) {
		std::fill(neighboursMap.begin(), neighboursMap.begin() + 9, -1);
	}
	if (_cluster[2] + 1 == _settings.layers) {
		std::fill(neighboursMap.begin() + 18, neighboursMap.begin() + 27, -1);
	}

	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	std::vector<int> sortedMap(neighboursMap);
	std::sort(sortedMap.begin(), sortedMap.end());

	for (size_t i = 0; i < sortedMap.size(); i++) {
		if (0 <= sortedMap[i] && sortedMap[i] < _settings.size) {
			eslocal index = std::find(neighboursMap.begin(), neighboursMap.end(), sortedMap[i]) - neighboursMap.begin();
			eslocal sx = index % 3 == 2 ? cNodes[0] - 1 : 0,
					ex = index % 3 != 0 ? cNodes[0] - 1 : 0,
					sy = index % 9 >= 6 ? cNodes[1] - 1 : 0,
					ey = index % 9 >= 3 ? cNodes[1] - 1 : 0,
					sz = index / 9 == 2 ? cNodes[2] - 1 : 0,
					ez = index / 9 >= 1 ? cNodes[2] - 1 : 0;

			for (eslocal x = sx; x <= ex; x++) {
				for (eslocal y = sy; y <= ey; y++) {
					for (eslocal z = sz; z <= ez; z++) {
						nodes[z * cNodes[0] * cNodes[1] + y * cNodes[0] + x]->clusters().push_back(sortedMap[i] + _settings.clusterOffset);
					}
				}
			}
			if (sortedMap[i] != config::env::MPIrank) {
				neighbours.push_back(sortedMap[i] + _settings.clusterOffset);
			}
		}
	}
}

}
}

