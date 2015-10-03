
#include "prisma15.h"

using namespace esinput;

eslocal Prisma15::subelements = Prisma15Subelements;

eslocal Prisma15::subnodes[3] = {
		Prisma15Subnodes,
		Prisma15Subnodes,
		Prisma15Subnodes,
};

Prisma15::Prisma15(const esinput::Settings &settings): _settings(settings)
{
	eslocal nodes[3];
	Utils<Prisma15>::clusterNodesCount(_settings, nodes);

	eslocal local = 0;
	for (eslocal z = 0; z < nodes[2]; z++) {
		for (eslocal y = 0; y < nodes[1]; y++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				_projection.push_back(local);
				if (addPoint(x, y, z)) {
					local++;
				}
			}
		}
	}

	Utils<Prisma6>::globalNodesCount(_settings, _g2Nodes);
	Utils<Prisma15>::globalNodesCount(_settings, _g3Nodes);
}

void Prisma15::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal prisma[20];
	prisma[0] = _projection[indices[0]];
	prisma[1] = _projection[indices[2]];
	prisma[2] = _projection[indices[8]];
	prisma[3] = _projection[indices[8]];
	prisma[4] = _projection[indices[18]];
	prisma[5] = _projection[indices[20]];
	prisma[6] = _projection[indices[26]];
	prisma[7] = _projection[indices[26]];

	prisma[8] = _projection[indices[1]];
	prisma[9] = _projection[indices[5]];
	prisma[10] = _projection[indices[8]];
	prisma[11] = _projection[indices[4]];
	prisma[12] = _projection[indices[19]];
	prisma[13] = _projection[indices[23]];
	prisma[14] = _projection[indices[26]];
	prisma[15] = _projection[indices[22]];
	prisma[16] = _projection[indices[9]];
	prisma[17] = _projection[indices[11]];
	prisma[18] = _projection[indices[17]];
	prisma[19] = _projection[indices[17]];
	mesh.pushElement(new mesh::Prisma15(prisma));

	prisma[0] = _projection[indices[0]];
	prisma[1] = _projection[indices[8]];
	prisma[2] = _projection[indices[6]];
	prisma[3] = _projection[indices[6]];
	prisma[4] = _projection[indices[18]];
	prisma[5] = _projection[indices[26]];
	prisma[6] = _projection[indices[24]];
	prisma[7] = _projection[indices[24]];

	prisma[8] = _projection[indices[4]];
	prisma[9] = _projection[indices[7]];
	prisma[10] = _projection[indices[6]];
	prisma[11] = _projection[indices[3]];
	prisma[12] = _projection[indices[22]];
	prisma[13] = _projection[indices[25]];
	prisma[14] = _projection[indices[24]];
	prisma[15] = _projection[indices[21]];
	prisma[16] = _projection[indices[9]];
	prisma[17] = _projection[indices[17]];
	prisma[18] = _projection[indices[15]];
	prisma[19] = _projection[indices[15]];
	mesh.pushElement(new mesh::Prisma15(prisma));
}

eslocal Prisma15::clusterNodesCount(const esinput::Settings &settings)
{
	eslocal nodes[3];
	Utils<Prisma15>::clusterNodesCount(settings, nodes);
	eslocal cElems = Utils<Prisma15>::clusterElementsCount(settings);

	eslocal count = nodes[0] * nodes[1] * nodes[2];	// Full mesh
	count -= 3 * (cElems / 2);						// remove unused nodes in mesh
	for (int i = 1; i < 3; i++) {					// remove unused nodes from the surface
		count -=
			settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] *
			settings.subdomainsInCluster[(i + 1) % 3] * settings.elementsInSubdomain[(i + 1) % 3];
	}

	return count;
}

esglobal Prisma15::globalNodesCount(const esinput::Settings &settings)
{
	esglobal nodes[3];
	Utils<Prisma15>::globalNodesCount(settings, nodes);
	esglobal cElems = Utils<Prisma15>::clusterElementsCount(settings);
	cElems *= settings.clusters[0] * settings.clusters[1] * settings.clusters[2];

	esglobal count = nodes[0] * nodes[1] * nodes[2];	// Full mesh
	count -= 3 * (cElems / 2);							// remove unused nodes in mesh
	for (int i = 1; i < 3; i++) {						// remove unused nodes from the surface
		count -=
			settings.clusters[i] * settings.clusters[(i + 1) % 3] *
			settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] *
			settings.subdomainsInCluster[(i + 1) % 3] * settings.elementsInSubdomain[(i + 1) % 3];
	}

	return count;
}



