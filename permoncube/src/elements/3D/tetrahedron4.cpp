#include "tetrahedron4.h"

using namespace permoncube;

size_t Tetrahedron4::subelements = Tetrahedron4Subelements;

size_t Tetrahedron4::subnodes[3] = {
		Tetrahedron4Subnodes,
		Tetrahedron4Subnodes,
		Tetrahedron4Subnodes
};

std::vector<idx_t> Tetrahedron4::_coordinateMapping;

void Tetrahedron4::addElements(mesh::Mesh &mesh, const idx_t indices[])
{
	idx_t tetra[5];
	tetra[0] = _coordinateMapping[indices[0]];
	tetra[1] = _coordinateMapping[indices[3]];
	tetra[2] = _coordinateMapping[indices[2]];
	tetra[4] = _coordinateMapping[indices[4]];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = _coordinateMapping[indices[3]];
	tetra[1] = _coordinateMapping[indices[2]];
	tetra[2] = _coordinateMapping[indices[4]];
	tetra[4] = _coordinateMapping[indices[6]];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = _coordinateMapping[indices[7]];
	tetra[1] = _coordinateMapping[indices[3]];
	tetra[2] = _coordinateMapping[indices[4]];
	tetra[4] = _coordinateMapping[indices[6]];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = _coordinateMapping[indices[3]];
	tetra[1] = _coordinateMapping[indices[5]];
	tetra[2] = _coordinateMapping[indices[7]];
	tetra[4] = _coordinateMapping[indices[4]];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = _coordinateMapping[indices[1]];
	tetra[1] = _coordinateMapping[indices[5]];
	tetra[2] = _coordinateMapping[indices[3]];
	tetra[4] = _coordinateMapping[indices[4]];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = _coordinateMapping[indices[0]];
	tetra[1] = _coordinateMapping[indices[4]];
	tetra[2] = _coordinateMapping[indices[1]];
	tetra[4] = _coordinateMapping[indices[3]];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));
}

void Tetrahedron4::addCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[])
{
	mesh::Coordinates &coordinates = mesh.coordinates();

	size_t nodes[3];

	ElementGenerator<Tetrahedron4>::clusterNodesCount(settings, nodes);
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

void Tetrahedron4::clear()
{
	_coordinateMapping.clear();
}

