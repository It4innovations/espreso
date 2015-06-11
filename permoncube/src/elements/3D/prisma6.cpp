
#include "prisma6.h"

using namespace permoncube;

eslocal Prisma6::subelements = Prisma6Subelements;

eslocal Prisma6::subnodes[3] = {
		Prisma6Subnodes,
		Prisma6Subnodes,
		Prisma6Subnodes,
};

void Prisma6::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal prisma[8];
	prisma[0] = indices[0];
	prisma[1] = indices[1];
	prisma[2] = indices[3];
	prisma[3] = indices[3];
	prisma[4] = indices[4];
	prisma[5] = indices[5];
	prisma[6] = indices[7];
	prisma[7] = indices[7];
	mesh.pushElement(new mesh::Prisma6(prisma));

	prisma[0] = indices[0];
	prisma[1] = indices[3];
	prisma[2] = indices[2];
	prisma[3] = indices[2];
	prisma[4] = indices[4];
	prisma[5] = indices[7];
	prisma[6] = indices[6];
	prisma[7] = indices[6];
	mesh.pushElement(new mesh::Prisma6(prisma));
}

void Prisma6::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	Element3D<Prisma6>::addFullCoordinates(mesh, settings, cluster);
}

void Prisma6::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Prisma6>::fixFullZeroPlanes(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Prisma6::fixBottom(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	Element3D<Prisma6>::fixFullBottom(settings, dirichlet_x, dirichlet_y, dirichlet_z, cluster);
}

void Prisma6::fillGlobalBoundaries(
		const permoncube::Settings &settings,
		mesh::Boundaries &boundaries)
{
	Element3D<Prisma6>::fillGlobalBoundaries(settings, boundaries);
}




