
#include "hexahedron20.h"

using namespace permoncube;

eslocal Hexahedron20::subelements = Hexahedron20Subelements;

eslocal Hexahedron20::subnodes[3] = {
		Hexahedron20Subnodes,
		Hexahedron20Subnodes,
		Hexahedron20Subnodes,
};

Hexahedron20::Hexahedron20(const permoncube::Settings &settings): _settings(settings)
{
	eslocal nodes[3];
	Utils<Hexahedron20>::clusterNodesCount(_settings, nodes);

	_projection.reserve(clusterNodesCount(_settings));
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
}

void Hexahedron20::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal hexa[20];
	hexa[0] = _projection[indices[2]];
	hexa[1] = _projection[indices[8]];
	hexa[2] = _projection[indices[6]];
	hexa[3] = _projection[indices[0]];
	hexa[4] = _projection[indices[20]];
	hexa[5] = _projection[indices[26]];
	hexa[6] = _projection[indices[24]];
	hexa[7] = _projection[indices[18]];

	hexa[8] = _projection[indices[5]];
	hexa[9] = _projection[indices[7]];
	hexa[10] = _projection[indices[3]];
	hexa[11] = _projection[indices[1]];
	hexa[12] = _projection[indices[23]];
	hexa[13] = _projection[indices[25]];
	hexa[14] = _projection[indices[21]];
	hexa[15] = _projection[indices[19]];
	hexa[16] = _projection[indices[11]];
	hexa[17] = _projection[indices[17]];
	hexa[18] = _projection[indices[15]];
	hexa[19] = _projection[indices[9]];
	mesh.pushElement(new mesh::Hexahedron20(hexa));
}

eslocal Hexahedron20::clusterNodesCount(const permoncube::Settings &settings)
{
	eslocal nodes[3];
	Utils<Hexahedron20>::clusterNodesCount(settings, nodes);
	eslocal cElems = Utils<Hexahedron20>::clusterElementsCount(settings);

	eslocal count = nodes[0] * nodes[1] * nodes[2];	// Full mesh
	count -= 4 * cElems;							// remove unused nodes in mesh
	for (int i = 0; i < 3; i++) {					// remove unused nodes from the surface
		count -=
			settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] *
			settings.subdomainsInCluster[(i + 1) % 3] * settings.elementsInSubdomain[(i + 1) % 3];
	}

	return count;
}

esglobal Hexahedron20::globalNodesCount(const permoncube::Settings &settings)
{
	esglobal nodes[3];
	Utils<Hexahedron20>::globalNodesCount(settings, nodes);
	esglobal cElems = Utils<Hexahedron20>::clusterElementsCount(settings);
	cElems *= settings.clusters[0] * settings.clusters[1] * settings.clusters[2];

	esglobal count = nodes[0] * nodes[1] * nodes[2];	// Full mesh
	count -= 4 * cElems;								// remove unused nodes in mesh
	for (int i = 0; i < 3; i++) {						// remove unused nodes from the surface
		count -=
			settings.clusters[i] * settings.clusters[(i + 1) % 3] *
			settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] *
			settings.subdomainsInCluster[(i + 1) % 3] * settings.elementsInSubdomain[(i + 1) % 3];
	}

	return count;
}


