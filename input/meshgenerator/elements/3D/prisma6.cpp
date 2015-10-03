
#include "../../../meshgenerator/elements/3D/prisma6.h"

using namespace esinput;

eslocal Prisma6::subelements = Prisma6Subelements;

eslocal Prisma6::subnodes[3] = {
		Prisma6Subnodes,
		Prisma6Subnodes,
		Prisma6Subnodes,
};

Prisma6::Prisma6(const esinput::CubeSettings &settings): _settings(settings)
{
	Utils<Prisma6>::globalNodesCount(_settings, _gNodes);
}

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

eslocal Prisma6::clusterNodesCount(const esinput::CubeSettings &settings)
{
	eslocal nodes[3];
	Utils<Prisma6>::clusterNodesCount(settings, nodes);

	return nodes[0] * nodes[1] * nodes[2];
}

esglobal Prisma6::globalNodesCount(const esinput::CubeSettings &settings)
{
	esglobal nodes[3];
	Utils<Prisma6>::globalNodesCount(settings, nodes);

	return nodes[0] * nodes[1] * nodes[2];
}


