
#include "pyramid5.h"

using namespace espreso;

Pyramid5Generator::Pyramid5Generator()
{
	subelements = 6;
	subnodes[0] = 3;
	subnodes[1] = 3;
	subnodes[2] = 3;
	enodes = 5;
	code = Element::CODE::PYRAMID5;
}

void Pyramid5Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
	elements.push_back(indices[18]);
	elements.push_back(indices[20]);
	elements.push_back(indices[ 2]);
	elements.push_back(indices[ 0]);
	elements.push_back(indices[13]);

	elements.push_back(indices[20]);
	elements.push_back(indices[26]);
	elements.push_back(indices[ 8]);
	elements.push_back(indices[ 2]);
	elements.push_back(indices[13]);

	elements.push_back(indices[26]);
	elements.push_back(indices[24]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[ 8]);
	elements.push_back(indices[13]);

	elements.push_back(indices[24]);
	elements.push_back(indices[18]);
	elements.push_back(indices[ 0]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[13]);

	elements.push_back(indices[18]);
	elements.push_back(indices[24]);
	elements.push_back(indices[26]);
	elements.push_back(indices[20]);
	elements.push_back(indices[13]);

	elements.push_back(indices[ 2]);
	elements.push_back(indices[ 8]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[ 0]);
	elements.push_back(indices[13]);
}

void Pyramid5Generator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const
{
	return;
}

void Pyramid5Generator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
	switch (face) {
	case CubeFace::X_1:
		nodes.push_back(indices[ 2]);
		nodes.push_back(indices[ 8]);
		nodes.push_back(indices[26]);
		nodes.push_back(indices[20]);
		break;
	case CubeFace::Y_1:
		nodes.push_back(indices[ 8]);
		nodes.push_back(indices[ 6]);
		nodes.push_back(indices[24]);
		nodes.push_back(indices[26]);
		break;
	case CubeFace::X_0:
		nodes.push_back(indices[ 6]);
		nodes.push_back(indices[ 0]);
		nodes.push_back(indices[18]);
		nodes.push_back(indices[24]);
		break;
	case CubeFace::Y_0:
		nodes.push_back(indices[ 0]);
		nodes.push_back(indices[ 2]);
		nodes.push_back(indices[20]);
		nodes.push_back(indices[18]);
		break;
	case CubeFace::Z_0:
		nodes.push_back(indices[ 0]);
		nodes.push_back(indices[ 6]);
		nodes.push_back(indices[ 8]);
		nodes.push_back(indices[ 2]);
		break;
	case CubeFace::Z_1:
		nodes.push_back(indices[20]);
		nodes.push_back(indices[26]);
		nodes.push_back(indices[24]);
		nodes.push_back(indices[18]);
		break;
	default:
		break;
	}
}

void Pyramid5Generator::pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const
{
	return;
//	pushNodes(elements, indices, edge);
//	esize.push_back(2);
//	etype.push_back((int)Element::CODE::LINE2);
}

void Pyramid5Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
	pushNodes(elements, indices, face);
	esize.push_back(4);
	etype.push_back((int)Element::CODE::SQUARE4);
}



