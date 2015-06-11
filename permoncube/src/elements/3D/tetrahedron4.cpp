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

eslocal Tetrahedron4::clusterNodesCount(const permoncube::Settings &settings)
{
	eslocal nodes[3];
	Utils<Tetrahedron4>::clusterNodesCount(settings, nodes);

	return nodes[0] * nodes[1] * nodes[2];
}

esglobal Tetrahedron4::globalNodesCount(const permoncube::Settings &settings)
{
	esglobal nodes[3];
	Utils<Tetrahedron4>::globalNodesCount(settings, nodes);

	return nodes[0] * nodes[1] * nodes[2];
}

