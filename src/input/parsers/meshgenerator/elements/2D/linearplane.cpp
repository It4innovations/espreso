
#include "linearplane.h"

using namespace espreso;

LinearPlaneGenerator::LinearPlaneGenerator()
{
	subnodes[0] = 2;
	subnodes[1] = 2;
	subnodes[2] = 1;
}

void LinearPlaneGenerator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const
{
	switch (edge) {
	case CubeEdge::X_0_Z_0:
		nodes.push_back(indices[2]);
		nodes.push_back(indices[0]);
		break;
	case CubeEdge::X_1_Z_0:
		nodes.push_back(indices[1]);
		nodes.push_back(indices[3]);
		break;
	case CubeEdge::Y_0_Z_0:
		nodes.push_back(indices[0]);
		nodes.push_back(indices[1]);
		break;
	case CubeEdge::Y_1_Z_0:
		nodes.push_back(indices[3]);
		nodes.push_back(indices[2]);
		break;
	default:
		return;
	}
}

void LinearPlaneGenerator::pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const
{
	pushNodes(elements, indices, edge);
	esize.push_back(2);
	etype.push_back((int)Element::CODE::LINE2);
}

void LinearPlaneGenerator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
	nodes.push_back(indices[0]);
	nodes.push_back(indices[1]);
	nodes.push_back(indices[2]);
	nodes.push_back(indices[3]);
	return;
}

void LinearPlaneGenerator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
	return;
}


