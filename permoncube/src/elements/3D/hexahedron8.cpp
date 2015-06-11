
#include "hexahedron8.h"

using namespace permoncube;

eslocal Hexahedron8::subelements = Hexahedron8Subelements;

eslocal Hexahedron8::subnodes[3] = {
		Hexahedron8Subnodes,
		Hexahedron8Subnodes,
		Hexahedron8Subnodes,
};

void Hexahedron8::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal hexa[8];
	hexa[0] = indices[0];
	hexa[1] = indices[1];
	hexa[2] = indices[3];
	hexa[3] = indices[2];
	hexa[4] = indices[4];
	hexa[5] = indices[5];
	hexa[6] = indices[7];
	hexa[7] = indices[6];
	mesh.pushElement(new mesh::Hexahedron8(hexa));
}

void Hexahedron8::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	Element3D<Hexahedron8>::addFullCoordinates(mesh, settings, cluster);
}

void Hexahedron8::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Hexahedron8>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Hexahedron8::fixBottom(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Hexahedron8>::fixFullBottom(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Hexahedron8::fillGlobalBoundaries(
		const permoncube::Settings &settings,
		mesh::Boundaries &boundaries)
{
	Element3D<Hexahedron8>::fillGlobalBoundaries(settings, boundaries);
}

