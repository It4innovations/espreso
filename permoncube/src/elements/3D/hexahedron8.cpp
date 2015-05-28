
#include "hexahedron8.h"

using namespace permoncube;

size_t Hexahedron8::subelements = Hexahedron8Subelements;

size_t Hexahedron8::subnodes[3] = {
		Hexahedron8Subnodes,
		Hexahedron8Subnodes,
		Hexahedron8Subnodes,
};

std::vector<idx_t> Hexahedron8::_coordinateMapping;

void Hexahedron8::addElements(mesh::Mesh &mesh, const idx_t indices[])
{
	idx_t hexa[8];
	hexa[0] = _coordinateMapping[indices[0]];
	hexa[1] = _coordinateMapping[indices[1]];
	hexa[2] = _coordinateMapping[indices[3]];
	hexa[3] = _coordinateMapping[indices[2]];
	hexa[4] = _coordinateMapping[indices[4]];
	hexa[5] = _coordinateMapping[indices[5]];
	hexa[6] = _coordinateMapping[indices[7]];
	hexa[7] = _coordinateMapping[indices[6]];
	mesh.pushElement(new mesh::Hexahedron8(hexa));
}

void Hexahedron8::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	mesh::Coordinates &coordinates = mesh.coordinates();

	size_t nodes[3];

	ElementGenerator<Hexahedron8>::clusterNodesCount(settings, nodes);
	mesh.coordinates().resize(nodes[0] * nodes[1] * nodes[2]);

	idx_t global = 0;
	idx_t local = 0;
	idx_t s[3], e[3];
	double step[3];
	for (int i = 0; i < 3; i++) {
		s[i] = (settings.subdomainsInCluster[i] * (settings.elementsInSubdomain[i] + subnodes[i])) * cluster[i];
		e[i] = (settings.subdomainsInCluster[i] * (settings.elementsInSubdomain[i] + subnodes[i])) * (cluster[i] + 1);
	}
	for (int i = 0; i < 3; i++) {
		step[i] = settings.clusterLength[i] / (nodes[i] - 1);
	}

	ElementGenerator<Hexahedron8>::globalNodesCount(settings, nodes);
	_coordinateMapping.resize(nodes[0] * nodes[1] * nodes[2]);

	for (idx_t z = 0; z < nodes[2]; z++) {
		for (idx_t y = 0; y < nodes[1]; y++) {
			for (idx_t x = 0; x < nodes[0]; x++) {
				_coordinateMapping[global] = local;
				if (s[2] <= z && z <= e[2] && s[1] <= y && y <= e[1] && s[0] <= x && x <= e[0]) {
					coordinates.add(global, mesh::Point(x * step[0], y * step[1], z * step[2]));
					local++;
				}
				global++;
			}
		}
	}
}

void Hexahedron8::fixZeroPlanes(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	size_t nodes[3];
	ElementGenerator<Hexahedron8>::globalNodesCount(settings, nodes);
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

void Hexahedron8::fixBottom(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	size_t nodes[3];
	ElementGenerator<Hexahedron8>::globalNodesCount(settings, nodes);
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

void Hexahedron8::clear()
{
	_coordinateMapping.clear();
}

