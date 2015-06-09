
#include "hexahedron20.h"

using namespace permoncube;

eslocal Hexahedron20::subelements = Hexahedron8Subelements;

eslocal Hexahedron20::subnodes[3] = {
		Hexahedron20Subnodes,
		Hexahedron20Subnodes,
		Hexahedron20Subnodes,
};

void Hexahedron20::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal hexa[20];
	hexa[0] = indices[2];
	hexa[1] = indices[8];
	hexa[2] = indices[6];
	hexa[3] = indices[0];
	hexa[4] = indices[20];
	hexa[5] = indices[26];
	hexa[6] = indices[24];
	hexa[7] = indices[18];

	hexa[8] = indices[8];
	hexa[9] = indices[7];
	hexa[10] = indices[3];
	hexa[11] = indices[1];
	hexa[12] = indices[23];
	hexa[13] = indices[25];
	hexa[14] = indices[21];
	hexa[15] = indices[19];
	hexa[16] = indices[11];
	hexa[17] = indices[17];
	hexa[18] = indices[15];
	hexa[19] = indices[9];
	mesh.pushElement(new mesh::Hexahedron20(hexa));
}

void Hexahedron20::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	Element3D<Hexahedron20>::addFullCoordinates(mesh, settings, cluster);
}

void Hexahedron20::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Hexahedron20>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Hexahedron20::fixBottom(
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
	Utils<Hexahedron20>::clusterNodesCount(settings, nodes);

	eslocal index = -1;
	for (eslocal y = 0; y < nodes[1]; y++) {
		for (eslocal x = 0; x < nodes[0]; x++) {
			index++;
			if ((y % 2 == 0) && (x % 2 == 0)) {
				continue;
			}
			dirichlet_z[index] = 0;
			dirichlet_y[index] = 0;
			dirichlet_x[index] = 0;
		}
	}
}

void Hexahedron20::fillGlobalBoundaries(
		const permoncube::Settings &settings,
		mesh::Boundaries &boundaries)
{
	Element3D<Hexahedron20>::fillGlobalBoundaries(settings, boundaries);
}


