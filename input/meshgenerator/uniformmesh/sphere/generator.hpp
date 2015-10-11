
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

	mesh::Point point;
	eslocal index = 0;
	for (eslocal z = 0; z < cNodes[2]; z++) {
		for (eslocal y = 0; y < cNodes[1]; y++) {
			for (eslocal x = 0; x < cNodes[0]; x++) {
				point.z = 1;
				point.y = (corner[0] * (1 - y * step[1]) + corner[2] * y * step[1]).y;
				point.x = (corner[0] * (1 - x * step[0]) + corner[1] * x * step[0]).x;

				point.normalize();
				point *= _settings.innerRadius + layerSize * (this->_rank / 6) + z * step[2];

				switch (this->_rank % 6) {
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

