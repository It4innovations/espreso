
#include "element3D.h"

using namespace permoncube;

template<class TElement>
void Element3D<TElement>::addFullCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[])
{
	mesh::Coordinates &coordinates = mesh.coordinates();

	eslocal cNodes[3];
	Utils<TElement>::clusterNodesCount(settings, cNodes);
	mesh.coordinates().reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	esglobal global = 0;
	eslocal local = 0;
	esglobal s[3], e[3];
	double step[3];
	for (eslocal i = 0; i < 3; i++) {
		s[i] = (cNodes[i] - 1) * cluster[i];
		e[i] = (cNodes[i] - 1) * (cluster[i] + 1);
	}
	for (eslocal i = 0; i < 3; i++) {
		step[i] = settings.clusterLength[i] / (cNodes[i] - 1);
	}

	esglobal gNodes[3];
	Utils<TElement>::globalNodesCount(settings, gNodes);

	for (esglobal z = 0; z < gNodes[2]; z++) {
		for (esglobal y = 0; y < gNodes[1]; y++) {
			for (esglobal x = 0; x < gNodes[0]; x++) {
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
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	if (cluster[2] > 0) {
		return;
	}
	eslocal nodes[3];
	Utils<TElement>::clusterNodesCount(settings, nodes);

	eslocal index = 0;
	for (eslocal y = 0; y < nodes[1]; y++) {
		for (eslocal x = 0; x < nodes[0]; x++) {
			dirichlet_z[index] = 0;
			dirichlet_y[index] = 0;
			dirichlet_x[index] = 0;
			index++;
		}
	}
}


template<class TElement>
void Element3D<TElement>::fixFullZeroPlanes(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	eslocal nodes[3];
	Utils<TElement>::clusterNodesCount(settings, nodes);

	if (cluster[0] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal y = 0; y < nodes[1]; y++) {
				dirichlet_x[index] = 0;
				index += nodes[0];
			}
		}
	}
	if (cluster[1] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				dirichlet_y[index] = 0;
				index++;
			}
			index += nodes[1] * nodes[0];
		}
	}
	if (cluster[2] == 0) {
		eslocal index = 0;
		for (eslocal y = 0; y < nodes[1]; y++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				dirichlet_z[index] = 0;
				index++;
			}
		}
	}
}


