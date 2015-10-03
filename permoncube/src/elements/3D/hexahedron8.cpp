
#include "hexahedron8.h"

using namespace esinput;

eslocal Hexahedron8::subelements = Hexahedron8Subelements;

eslocal Hexahedron8::subnodes[3] = {
		Hexahedron8Subnodes,
		Hexahedron8Subnodes,
		Hexahedron8Subnodes,
};

Hexahedron8::Hexahedron8(const esinput::Settings &settings): _settings(settings)
{
	Utils<Hexahedron8>::globalNodesCount(_settings, _gNodes);
}

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

eslocal Hexahedron8::clusterNodesCount(const esinput::Settings &settings)
{
	eslocal nodes[3];
	Utils<Hexahedron8>::clusterNodesCount(settings, nodes);

	return nodes[0] * nodes[1] * nodes[2];
}

esglobal Hexahedron8::globalNodesCount(const esinput::Settings &settings)
{
	esglobal nodes[3];
	Utils<Hexahedron8>::globalNodesCount(settings, nodes);

	return nodes[0] * nodes[1] * nodes[2];
}

