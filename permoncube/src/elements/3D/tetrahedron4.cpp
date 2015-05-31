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
	Element3D<Tetrahedron4>::addFullCoordinates(mesh, settings, cluster, _coordinateMapping);
}

void Tetrahedron4::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<int, double> &dirichlet_x,
		std::map<int, double> &dirichlet_y,
		std::map<int, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Tetrahedron4>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster, _coordinateMapping);
}

void Tetrahedron4::fixBottom(
		const permoncube::Settings &settings,
		std::map<int, double> &dirichlet_x,
		std::map<int, double> &dirichlet_y,
		std::map<int, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Tetrahedron4>::fixFullBottom(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster, _coordinateMapping);
}

void Tetrahedron4::clear()
{
	_coordinateMapping.clear();
}

