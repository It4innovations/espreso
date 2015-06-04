
#include "tetrahedron10.h"

using namespace permoncube;

eslocal Tetrahedron10::subelements = Tetrahedron10Subelements;

eslocal Tetrahedron10::subnodes[3] = {
		Tetrahedron10Subnodes,
		Tetrahedron10Subnodes,
		Tetrahedron10Subnodes
};

void Tetrahedron10::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal tetra[20];
	tetra[0] = indices[2];
	tetra[1] = indices[6];
	tetra[2] = indices[0];
	tetra[4] = indices[20];

	tetra[8] = indices[4];
	tetra[9] = indices[3];
	tetra[11] = indices[1];
	tetra[16] = indices[11];
	tetra[17] = indices[13];
	tetra[18] = indices[10];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = indices[6];
	tetra[1] = indices[0];
	tetra[2] = indices[20];
	tetra[4] = indices[18];

	tetra[8] = indices[3];
	tetra[9] = indices[10];
	tetra[11] = indices[13];
	tetra[16] = indices[12];
	tetra[17] = indices[9];
	tetra[18] = indices[19];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = indices[24];
	tetra[1] = indices[6];
	tetra[2] = indices[20];
	tetra[4] = indices[18];

	tetra[8] = indices[15];
	tetra[9] = indices[13];
	tetra[11] = indices[22];
	tetra[16] = indices[21];
	tetra[17] = indices[12];
	tetra[18] = indices[19];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = indices[6];
	tetra[1] = indices[26];
	tetra[2] = indices[24];
	tetra[4] = indices[20];

	tetra[8] = indices[16];
	tetra[9] = indices[25];
	tetra[11] = indices[15];
	tetra[16] = indices[13];
	tetra[17] = indices[23];
	tetra[18] = indices[22];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = indices[8];
	tetra[1] = indices[26];
	tetra[2] = indices[6];
	tetra[4] = indices[20];

	tetra[8] = indices[17];
	tetra[9] = indices[16];
	tetra[11] = indices[7];
	tetra[16] = indices[14];
	tetra[17] = indices[23];
	tetra[18] = indices[13];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = indices[2];
	tetra[1] = indices[20];
	tetra[2] = indices[8];
	tetra[4] = indices[6];

	tetra[8] = indices[11];
	tetra[9] = indices[14];
	tetra[11] = indices[5];
	tetra[16] = indices[4];
	tetra[17] = indices[13];
	tetra[18] = indices[7];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));
}

void Tetrahedron10::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	Element3D<Tetrahedron10>::addFullCoordinates(mesh, settings, cluster);
}

void Tetrahedron10::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Tetrahedron10>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Tetrahedron10::fixBottom(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Tetrahedron10>::fixFullBottom(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Tetrahedron10::fillGlobalBoundaries(
		const permoncube::Settings &settings,
		mesh::Boundaries &boundaries)
{
	Element3D<Tetrahedron10>::fillGlobalBoundaries(settings, boundaries);
}



