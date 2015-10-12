
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
	: UniformGenerator<TElement>(argc, argv, rank, size), _settings(argc, argv)
{
	if (_settings.layers * 6 != this->_size) {
		if (this->_rank == 0) {
			std::cerr << "The number of clusters(" << _settings.layers * 6;
			std::cerr << ") does not accord the number of MPI processes(" << this->_size << ").\n";
		}
		exit(EXIT_FAILURE);
	}
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
				point *= _settings.innerRadius + layerSize * (this->_rank / 6) + z * step[2];

				switch (this->_rank % 6) {
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
	mesh::CoordinatesProperty &dirichlet_x = coordinates.property(mesh::DIRICHLET_X);
	mesh::CoordinatesProperty &dirichlet_y = coordinates.property(mesh::DIRICHLET_Y);
	mesh::CoordinatesProperty &dirichlet_z = coordinates.property(mesh::DIRICHLET_Z);

	eslocal surface = SphereUtils<TElement>::surfaceNodesCount(_settings);

	for (eslocal i = 0; i < surface; i++) {
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

	eslocal surface = SphereUtils<TElement>::surfaceNodesCount(_settings);
	if (_rank / 6 > 0) {
		for (eslocal i = 0; i < cNodes[1] * cNodes[0]; i++) {
			boundaries[i].insert(_rank - 6);
		}
	}

	if (_rank / 6 - 1 < _settings.layers) {
		for (eslocal i = 0; i < cNodes[1] * cNodes[0]; i++) {
			boundaries[i + (cNodes[2] - 1) * cNodes[1] * cNodes[0]].insert(_rank + 6);
		}
	}

	eslocal index;
	if  (_rank % 6 < 4) {
		index = 0;
		for (eslocal z = 0; z < cNodes[2]; z++) {
			for (eslocal y = 0; y < cNodes[1]; y++) {
				boundaries[index++].insert()
			}
		}
	}
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

