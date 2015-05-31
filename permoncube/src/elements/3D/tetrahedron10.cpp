
#include "tetrahedron10.h"

using namespace permoncube;

size_t Tetrahedron10::subelements = Tetrahedron10Subelements;

size_t Tetrahedron10::subnodes[3] = {
		Tetrahedron10Subnodes,
		Tetrahedron10Subnodes,
		Tetrahedron10Subnodes
};

std::vector<idx_t> Tetrahedron10::_coordinateMapping;

void Tetrahedron10::addElements(mesh::Mesh &mesh, const idx_t indices[])
{
	idx_t tetra[20];
	tetra[0] = _coordinateMapping[indices[2]];
	tetra[1] = _coordinateMapping[indices[6]];
	tetra[2] = _coordinateMapping[indices[0]];
	tetra[4] = _coordinateMapping[indices[20]];

	tetra[8] = _coordinateMapping[indices[4]];
	tetra[9] = _coordinateMapping[indices[3]];
	tetra[11] = _coordinateMapping[indices[1]];
	tetra[16] = _coordinateMapping[indices[11]];
	tetra[17] = _coordinateMapping[indices[13]];
	tetra[18] = _coordinateMapping[indices[10]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[6]];
	tetra[1] = _coordinateMapping[indices[0]];
	tetra[2] = _coordinateMapping[indices[20]];
	tetra[4] = _coordinateMapping[indices[18]];

	tetra[8] = _coordinateMapping[indices[3]];
	tetra[9] = _coordinateMapping[indices[10]];
	tetra[11] = _coordinateMapping[indices[13]];
	tetra[16] = _coordinateMapping[indices[12]];
	tetra[17] = _coordinateMapping[indices[9]];
	tetra[18] = _coordinateMapping[indices[19]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[24]];
	tetra[1] = _coordinateMapping[indices[6]];
	tetra[2] = _coordinateMapping[indices[20]];
	tetra[4] = _coordinateMapping[indices[18]];

	tetra[8] = _coordinateMapping[indices[15]];
	tetra[9] = _coordinateMapping[indices[13]];
	tetra[11] = _coordinateMapping[indices[22]];
	tetra[16] = _coordinateMapping[indices[21]];
	tetra[17] = _coordinateMapping[indices[12]];
	tetra[18] = _coordinateMapping[indices[19]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[6]];
	tetra[1] = _coordinateMapping[indices[26]];
	tetra[2] = _coordinateMapping[indices[24]];
	tetra[4] = _coordinateMapping[indices[20]];

	tetra[8] = _coordinateMapping[indices[16]];
	tetra[9] = _coordinateMapping[indices[25]];
	tetra[11] = _coordinateMapping[indices[15]];
	tetra[16] = _coordinateMapping[indices[13]];
	tetra[17] = _coordinateMapping[indices[23]];
	tetra[18] = _coordinateMapping[indices[22]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[8]];
	tetra[1] = _coordinateMapping[indices[26]];
	tetra[2] = _coordinateMapping[indices[6]];
	tetra[4] = _coordinateMapping[indices[20]];

	tetra[8] = _coordinateMapping[indices[17]];
	tetra[9] = _coordinateMapping[indices[16]];
	tetra[11] = _coordinateMapping[indices[7]];
	tetra[16] = _coordinateMapping[indices[14]];
	tetra[17] = _coordinateMapping[indices[23]];
	tetra[18] = _coordinateMapping[indices[13]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[2]];
	tetra[1] = _coordinateMapping[indices[20]];
	tetra[2] = _coordinateMapping[indices[8]];
	tetra[4] = _coordinateMapping[indices[6]];

	tetra[8] = _coordinateMapping[indices[11]];
	tetra[9] = _coordinateMapping[indices[14]];
	tetra[11] = _coordinateMapping[indices[5]];
	tetra[16] = _coordinateMapping[indices[4]];
	tetra[17] = _coordinateMapping[indices[13]];
	tetra[18] = _coordinateMapping[indices[7]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));
}

void Tetrahedron10::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	Element3D<Tetrahedron10>::addFullCoordinates(mesh, settings, cluster, _coordinateMapping);
}

void Tetrahedron10::fixZeroPlanes(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	Element3D<Tetrahedron10>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z);
}

void Tetrahedron10::fixBottom(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	Element3D<Tetrahedron10>::fixFullBottom(settings, dirichlet_x, dirichlet_y, dirichlet_z);
}

void Tetrahedron10::clear()
{
	_coordinateMapping.clear();
}


