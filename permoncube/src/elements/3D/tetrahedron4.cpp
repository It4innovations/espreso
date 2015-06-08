#include "tetrahedron4.h"

using namespace permoncube;

eslocal Tetrahedron4::subelements = Tetrahedron4Subelements;

eslocal Tetrahedron4::subnodes[3] = {
		Tetrahedron4Subnodes,
		Tetrahedron4Subnodes,
		Tetrahedron4Subnodes
};

void Tetrahedron4::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal tetra[5];
	tetra[0] = indices[0];
	tetra[1] = indices[3];
	tetra[2] = indices[2];
	tetra[4] = indices[4];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[3];
	tetra[1] = indices[2];
	tetra[2] = indices[4];
	tetra[4] = indices[6];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[7];
	tetra[1] = indices[3];
	tetra[2] = indices[4];
	tetra[4] = indices[6];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[3];
	tetra[1] = indices[5];
	tetra[2] = indices[7];
	tetra[4] = indices[4];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[1];
	tetra[1] = indices[5];
	tetra[2] = indices[3];
	tetra[4] = indices[4];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[0];
	tetra[1] = indices[4];
	tetra[2] = indices[1];
	tetra[4] = indices[3];
	mesh.pushElement(new mesh::Tetrahedron4(tetra));
}

void Tetrahedron4::addCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[])
{
	Element3D<Tetrahedron4>::addFullCoordinates(mesh, settings, cluster);
}

void Tetrahedron4::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Tetrahedron4>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Tetrahedron4::fixBottom(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Tetrahedron4>::fixFullBottom(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Tetrahedron4::fillGlobalBoundaries(
		const permoncube::Settings &settings,
		mesh::Boundaries &boundaries)
{
	Element3D<Tetrahedron4>::fillGlobalBoundaries(settings, boundaries);
}


