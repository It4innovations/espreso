
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
	Element3D<Tetrahedron4>::addFullCoordinates(mesh, settings, cluster, _coordinateMapping);
}

void Hexahedron8::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<int, double> &dirichlet_x,
		std::map<int, double> &dirichlet_y,
		std::map<int, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Hexahedron8>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster, _coordinateMapping);
}

void Hexahedron8::fixBottom(
		const permoncube::Settings &settings,
		std::map<int, double> &dirichlet_x,
		std::map<int, double> &dirichlet_y,
		std::map<int, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Hexahedron8>::fixFullBottom(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster, _coordinateMapping);
}

void Hexahedron8::clear()
{
	_coordinateMapping.clear();
}

