
#include "pyramid5.h"

using namespace permoncube;

eslocal Pyramid5::subelements = Pyramid5Subelements;

eslocal Pyramid5::subnodes[3] = {
		Pyramid5Subnodes,
		Pyramid5Subnodes,
		Pyramid5Subnodes,
};

Pyramid5::Pyramid5(const permoncube::Settings &settings): _settings(settings)
{
	eslocal nodes[3];
	Utils<Pyramid5>::clusterNodesCount(_settings, nodes);

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

void Pyramid5::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal pyramid[8];
	pyramid[0] = _projection[indices[18]];
	pyramid[1] = _projection[indices[20]];
	pyramid[2] = _projection[indices[2]];
	pyramid[3] = _projection[indices[0]];
	pyramid[4] = _projection[indices[13]];
	pyramid[5] = _projection[indices[13]];
	pyramid[6] = _projection[indices[13]];
	pyramid[7] = _projection[indices[13]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[20]];
	pyramid[1] = _projection[indices[26]];
	pyramid[2] = _projection[indices[8]];
	pyramid[3] = _projection[indices[2]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[26]];
	pyramid[1] = _projection[indices[24]];
	pyramid[2] = _projection[indices[6]];
	pyramid[3] = _projection[indices[8]];;
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[24]];
	pyramid[1] = _projection[indices[18]];
	pyramid[2] = _projection[indices[0]];
	pyramid[3] = _projection[indices[6]];;
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[18]];
	pyramid[1] = _projection[indices[24]];
	pyramid[2] = _projection[indices[26]];
	pyramid[3] = _projection[indices[20]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[2]];
	pyramid[1] = _projection[indices[8]];
	pyramid[2] = _projection[indices[6]];
	pyramid[3] = _projection[indices[0]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));
}

eslocal Pyramid5::clusterNodesCount(const permoncube::Settings &settings)
{
	Pyramid5::subelements = 1;
	for (int i = 0; i < 3; i++) {
		Pyramid5::subnodes[0] = 0;
	}

	eslocal nodes[3];
	Utils<Pyramid5>::clusterNodesCount(settings, nodes);
	eslocal cElems = Utils<Pyramid5>::clusterElementsCount(settings);

	eslocal count = nodes[0] * nodes[1] * nodes[2] + cElems;

	Pyramid5::subelements = Pyramid5Subelements;
	for (int i = 0; i < 3; i++) {
		Pyramid5::subnodes[0] = Pyramid5Subnodes;
	}

	return count;
}

esglobal Pyramid5::globalNodesCount(const permoncube::Settings &settings)
{
	Pyramid5::subelements = 1;
	for (int i = 0; i < 3; i++) {
		Pyramid5::subnodes[0] = 0;
	}

	esglobal nodes[3];
	Utils<Pyramid5>::globalNodesCount(settings, nodes);
	esglobal cElems = Utils<Pyramid5>::clusterElementsCount(settings);
	cElems *= settings.clusters[0] * settings.clusters[1] * settings.clusters[2];

	esglobal count = nodes[0] * nodes[1] * nodes[2] + cElems;	// Full mesh

	Pyramid5::subelements = Pyramid5Subelements;
	for (int i = 0; i < 3; i++) {
		Pyramid5::subnodes[0] = Pyramid5Subnodes;
	}

	return count;
}



