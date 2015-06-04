
#include "element3D.h"

using namespace permoncube;

template<class TElement>
void Element3D<TElement>::addFullCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[])
{
	mesh::Coordinates &coordinates = mesh.coordinates();

	esint cNodes[3];
	Utils<TElement>::clusterNodesCount(settings, cNodes);
	mesh.coordinates().reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	eslong global = 0;
	esint local = 0;
	eslong s[3], e[3];
	double step[3];
	for (esint i = 0; i < 3; i++) {
		s[i] = (cNodes[i] - 1) * cluster[i];
		e[i] = (cNodes[i] - 1) * (cluster[i] + 1);
	}
	for (esint i = 0; i < 3; i++) {
		step[i] = settings.clusterLength[i] / (cNodes[i] - 1);
	}

	eslong gNodes[3];
	Utils<TElement>::globalNodesCount(settings, gNodes);

	for (eslong z = 0; z < gNodes[2]; z++) {
		for (eslong y = 0; y < gNodes[1]; y++) {
			for (eslong x = 0; x < gNodes[0]; x++) {
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
		std::map<esint, double> &dirichlet_x,
		std::map<esint, double> &dirichlet_y,
		std::map<esint, double> &dirichlet_z,
		const size_t cluster[])
{
	if (cluster[2] > 0) {
		return;
	}
	esint nodes[3];
	Utils<TElement>::clusterNodesCount(settings, nodes);

	esint index = 0;
	for (esint y = 0; y < nodes[1]; y++) {
		for (esint x = 0; x < nodes[0]; x++) {
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
		std::map<esint, double> &dirichlet_x,
		std::map<esint, double> &dirichlet_y,
		std::map<esint, double> &dirichlet_z,
		const size_t cluster[])
{
	esint nodes[3];
	Utils<TElement>::clusterNodesCount(settings, nodes);

	if (cluster[0] == 0) {
		esint index = 0;
		for (esint z = 0; z < nodes[2]; z++) {
			for (esint y = 0; y < nodes[1]; y++) {
				dirichlet_x[index] = 0;
				index += nodes[0];
			}
		}
	}
	if (cluster[1] == 0) {
		esint index = 0;
		for (esint z = 0; z < nodes[2]; z++) {
			for (esint x = 0; x < nodes[0]; x++) {
				dirichlet_y[index] = 0;
				index++;
			}
			index += nodes[1] * nodes[0];
		}
	}
	if (cluster[2] == 0) {
		esint index = 0;
		for (esint y = 0; y < nodes[1]; y++) {
			for (esint x = 0; x < nodes[0]; x++) {
				dirichlet_z[index] = 0;
				index++;
			}
		}
	}
}


