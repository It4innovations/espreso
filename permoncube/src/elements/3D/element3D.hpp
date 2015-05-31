
#include "element3D.h"

using namespace permoncube;

template<class TElement>
void Element3D<TElement>::addFullCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[], std::vector<idx_t> &mapping)
{
	mesh::Coordinates &coordinates = mesh.coordinates();

	size_t nodes[3];

	Utils<TElement>::clusterNodesCount(settings, nodes);
	mesh.coordinates().reserve(nodes[0] * nodes[1] * nodes[2]);

	idx_t global = 0;
	idx_t local = 0;
	idx_t s[3], e[3];
	double step[3];
	for (int i = 0; i < 3; i++) {
		s[i] = (nodes[i] - 1) * cluster[i];
		e[i] = (nodes[i] - 1) * (cluster[i] + 1);
	}
	for (int i = 0; i < 3; i++) {
		step[i] = settings.clusterLength[i] / (nodes[i] - 1);
	}

	Utils<TElement>::globalNodesCount(settings, nodes);
	mapping.resize(nodes[0] * nodes[1] * nodes[2]);

	for (idx_t z = 0; z < nodes[2]; z++) {
		for (idx_t y = 0; y < nodes[1]; y++) {
			for (idx_t x = 0; x < nodes[0]; x++) {
				mapping[global] = local;
				if (s[2] <= z && z <= e[2] && s[1] <= y && y <= e[1] && s[0] <= x && x <= e[0]) {
					coordinates.add(mesh::Point(x * step[0], y * step[1], z * step[2]), local, global);
					local++;
				}
				global++;
			}
		}
	}
}

template<class TElement>
void Element3D<TElement>::fixFullBottom(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	size_t nodes[3];
	Utils<TElement>::globalNodesCount(settings, nodes);
	idx_t index = 0;
	for (idx_t y = 0; y < nodes[1]; y++) {
		for (idx_t x = 0; x < nodes[0]; x++) {
			dirichlet_z[index] = 0.0;
			dirichlet_y[index] = 0.0;
			dirichlet_x[index] = 0.0;
			index++;
		}
	}
}


template<class TElement>
void Element3D<TElement>::fixFullZeroPlanes(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	size_t nodes[3];
	Utils<TElement>::globalNodesCount(settings, nodes);
	idx_t index = 0;
	for (idx_t z = 0; z < nodes[2]; z++) {
		for (idx_t y = 0; y < nodes[1]; y++) {
			for (idx_t x = 0; x < nodes[0]; x++) {
				if (z == 0) {
					dirichlet_z[index] = 0.0;
				}
				if (y == 0) {
					dirichlet_y[index] = 0.0;
				}
				if (x == 0) {
					dirichlet_x[index] = 0.0;
				}
				index++;
			}
		}
	}
}


