
#include "../../../meshgenerator/elements/3D/pyramid13.h"

using namespace esinput;

eslocal Pyramid13::subelements = Pyramid13Subelements;

eslocal Pyramid13::subnodes[3] = {
		Pyramid13Subnodes,
		Pyramid13Subnodes,
		Pyramid13Subnodes,
};

Pyramid13::Pyramid13(const esinput::CubeSettings &settings): _settings(settings)
{
	eslocal nodes[3];
	Utils<Pyramid13>::clusterNodesCount(_settings, nodes);

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

	// Use elements with less sub-node -> it simplifies offset methods
	Utils<Hexahedron20>::globalNodesCount(_settings, _g3Nodes);
	Utils<Hexahedron8>::globalNodesCount(_settings, _g2Nodes);
	faceNodes = _g3Nodes[0] * _g3Nodes[1]  - (_g2Nodes[0] - 1) * (_g2Nodes[1] - 1);
}

void Pyramid13::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal pyramid[20];
	pyramid[0] = _projection[indices[100]];
	pyramid[1] = _projection[indices[104]];
	pyramid[2] = _projection[indices[4]];
	pyramid[3] = _projection[indices[0]];
	pyramid[4] = _projection[indices[62]];
	pyramid[5] = _projection[indices[62]];
	pyramid[6] = _projection[indices[62]];
	pyramid[7] = _projection[indices[62]];

	pyramid[8] = _projection[indices[102]];
	pyramid[9] = _projection[indices[54]];
	pyramid[10] = _projection[indices[2]];
	pyramid[11] = _projection[indices[50]];
	pyramid[12] = _projection[indices[62]];
	pyramid[13] = _projection[indices[62]];
	pyramid[14] = _projection[indices[62]];
	pyramid[15] = _projection[indices[62]];
	pyramid[16] = _projection[indices[81]];
	pyramid[17] = _projection[indices[83]];
	pyramid[18] = _projection[indices[33]];
	pyramid[19] = _projection[indices[31]];
	mesh.pushElement(new mesh::Pyramid13(pyramid));

	pyramid[0] = _projection[indices[104]];
	pyramid[1] = _projection[indices[124]];
	pyramid[2] = _projection[indices[24]];
	pyramid[3] = _projection[indices[4]];

	pyramid[8] = _projection[indices[114]];
	pyramid[9] = _projection[indices[74]];
	pyramid[10] = _projection[indices[14]];
	pyramid[11] = _projection[indices[54]];
	pyramid[16] = _projection[indices[83]];
	pyramid[17] = _projection[indices[93]];
	pyramid[18] = _projection[indices[43]];
	pyramid[19] = _projection[indices[33]];
	mesh.pushElement(new mesh::Pyramid13(pyramid));

	pyramid[0] = _projection[indices[124]];
	pyramid[1] = _projection[indices[120]];
	pyramid[2] = _projection[indices[20]];
	pyramid[3] = _projection[indices[24]];

	pyramid[8] = _projection[indices[122]];
	pyramid[9] = _projection[indices[70]];
	pyramid[10] = _projection[indices[22]];
	pyramid[11] = _projection[indices[74]];
	pyramid[16] = _projection[indices[93]];
	pyramid[17] = _projection[indices[91]];
	pyramid[18] = _projection[indices[41]];
	pyramid[19] = _projection[indices[43]];
	mesh.pushElement(new mesh::Pyramid13(pyramid));

	pyramid[0] = _projection[indices[120]];
	pyramid[1] = _projection[indices[100]];
	pyramid[2] = _projection[indices[0]];
	pyramid[3] = _projection[indices[20]];

	pyramid[8] = _projection[indices[110]];
	pyramid[9] = _projection[indices[50]];
	pyramid[10] = _projection[indices[10]];
	pyramid[11] = _projection[indices[70]];
	pyramid[16] = _projection[indices[91]];
	pyramid[17] = _projection[indices[81]];
	pyramid[18] = _projection[indices[31]];
	pyramid[19] = _projection[indices[41]];
	mesh.pushElement(new mesh::Pyramid13(pyramid));

	pyramid[0] = _projection[indices[100]];
	pyramid[1] = _projection[indices[120]];
	pyramid[2] = _projection[indices[124]];
	pyramid[3] = _projection[indices[104]];

	pyramid[8] = _projection[indices[110]];
	pyramid[9] = _projection[indices[122]];
	pyramid[10] = _projection[indices[114]];
	pyramid[11] = _projection[indices[102]];
	pyramid[16] = _projection[indices[81]];
	pyramid[17] = _projection[indices[91]];
	pyramid[18] = _projection[indices[93]];
	pyramid[19] = _projection[indices[83]];
	mesh.pushElement(new mesh::Pyramid13(pyramid));

	pyramid[0] = _projection[indices[4]];
	pyramid[1] = _projection[indices[24]];
	pyramid[2] = _projection[indices[20]];
	pyramid[3] = _projection[indices[0]];

	pyramid[8] = _projection[indices[14]];
	pyramid[9] = _projection[indices[22]];
	pyramid[10] = _projection[indices[10]];
	pyramid[11] = _projection[indices[2]];
	pyramid[16] = _projection[indices[33]];
	pyramid[17] = _projection[indices[43]];
	pyramid[18] = _projection[indices[41]];
	pyramid[19] = _projection[indices[31]];
	mesh.pushElement(new mesh::Pyramid13(pyramid));
}

eslocal Pyramid13::clusterNodesCount(const esinput::CubeSettings &settings)
{
	eslocal count = Hexahedron20::clusterNodesCount(settings);
	eslocal cElems = Utils<Hexahedron20>::clusterElementsCount(settings);

	return count + 9 * cElems;
}

esglobal Pyramid13::globalNodesCount(const esinput::CubeSettings &settings)
{
	eslocal count = Hexahedron20::globalNodesCount(settings);
	eslocal cElems = Utils<Hexahedron20>::clusterElementsCount(settings);
	cElems *= settings.clusters[0] * settings.clusters[1] * settings.clusters[2];

	return count + 0 * cElems;
}



