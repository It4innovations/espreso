
#include "linearvolume.h"

using namespace espreso;

LinearVolumeGenerator::LinearVolumeGenerator()
{
	subnodes[0] = 2;
	subnodes[1] = 2;
	subnodes[2] = 2;
}

void LinearVolumeGenerator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const
{
//	switch (edge) {
//	case CubeEdge::X_0_Y_0:
//		nodes.push_back(indices[0]);
//		nodes.push_back(indices[4]);
//		break;
//	case CubeEdge::X_0_Y_1:
//		nodes.push_back(indices[2]);
//		nodes.push_back(indices[6]);
//		break;
//	case CubeEdge::X_0_Z_0:
//		nodes.push_back(indices[2]);
//		nodes.push_back(indices[6]);
//		break;
//	case CubeEdge::X_0_Z_1:
//		nodes.push_back(indices[0]);
//		nodes.push_back(indices[4]);
//		break;
//	case CubeEdge::X_1_Y_0:
//		nodes.push_back(indices[0]);
//		nodes.push_back(indices[1]);
//		break;
//	case CubeEdge::X_1_Y_1:
//		nodes.push_back(indices[5]);
//		nodes.push_back(indices[4]);
//		break;
//	case CubeEdge::X_1_Z_0:
//		nodes.push_back(indices[1]);
//		nodes.push_back(indices[3]);
//		break;
//	case CubeEdge::X_1_Z_1:
//		nodes.push_back(indices[3]);
//		nodes.push_back(indices[7]);
//		break;
//	case CubeEdge::Y_0_Z_0:
//		nodes.push_back(indices[2]);
//		nodes.push_back(indices[6]);
//		break;
//	case CubeEdge::Y_0_Z_1:
//		nodes.push_back(indices[0]);
//		nodes.push_back(indices[4]);
//		break;
//	case CubeEdge::Y_1_Z_0:
//		nodes.push_back(indices[0]);
//		nodes.push_back(indices[1]);
//		break;
//	case CubeEdge::Y_1_Z_1:
//		nodes.push_back(indices[5]);
//		nodes.push_back(indices[4]);
//		break;
//	default:
//		return;
//	}
}

void LinearVolumeGenerator::pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const
{
	return;
//	pushNodes(elements, indices, edge);
//	esize.push_back(2);
//	etype.push_back((int)Element::CODE::LINE2);
}

void LinearVolumeGenerator::pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{

}

void LinearVolumeGenerator::pushSquareNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
	switch (face) {
	case CubeFace::X_1:
		nodes.push_back(indices[1]);
		nodes.push_back(indices[3]);
		nodes.push_back(indices[7]);
		nodes.push_back(indices[5]);
		break;
	case CubeFace::Y_1:
		nodes.push_back(indices[3]);
		nodes.push_back(indices[2]);
		nodes.push_back(indices[6]);
		nodes.push_back(indices[7]);
		break;
	case CubeFace::X_0:
		nodes.push_back(indices[2]);
		nodes.push_back(indices[0]);
		nodes.push_back(indices[4]);
		nodes.push_back(indices[6]);
		break;
	case CubeFace::Y_0:
		nodes.push_back(indices[0]);
		nodes.push_back(indices[1]);
		nodes.push_back(indices[5]);
		nodes.push_back(indices[4]);
		break;
	case CubeFace::Z_0:
		nodes.push_back(indices[0]);
		nodes.push_back(indices[2]);
		nodes.push_back(indices[3]);
		nodes.push_back(indices[1]);
		break;
	case CubeFace::Z_1:
		nodes.push_back(indices[5]);
		nodes.push_back(indices[7]);
		nodes.push_back(indices[6]);
		nodes.push_back(indices[4]);
		break;
	default:
		return;
	}
}

void LinearVolumeGenerator::pushTriangleNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
	switch (face) {
	case CubeFace::X_0:
		nodes.push_back(indices[2]);
		nodes.push_back(indices[0]);
		nodes.push_back(indices[4]);

		nodes.push_back(indices[2]);
		nodes.push_back(indices[4]);
		nodes.push_back(indices[6]);
		break;
	case CubeFace::X_1:
		nodes.push_back(indices[1]);
		nodes.push_back(indices[3]);
		nodes.push_back(indices[5]);

		nodes.push_back(indices[3]);
		nodes.push_back(indices[7]);
		nodes.push_back(indices[5]);
		break;
	case CubeFace::Y_0:
		nodes.push_back(indices[0]);
		nodes.push_back(indices[1]);
		nodes.push_back(indices[4]);

		nodes.push_back(indices[1]);
		nodes.push_back(indices[5]);
		nodes.push_back(indices[4]);
		break;
	case CubeFace::Y_1:
		nodes.push_back(indices[3]);
		nodes.push_back(indices[2]);
		nodes.push_back(indices[6]);

		nodes.push_back(indices[3]);
		nodes.push_back(indices[6]);
		nodes.push_back(indices[7]);
		break;
	case CubeFace::Z_0:
		nodes.push_back(indices[0]);
		nodes.push_back(indices[2]);
		nodes.push_back(indices[3]);

		nodes.push_back(indices[0]);
		nodes.push_back(indices[3]);
		nodes.push_back(indices[1]);
		break;
	case CubeFace::Z_1:
		nodes.push_back(indices[5]);
		nodes.push_back(indices[7]);
		nodes.push_back(indices[4]);

		nodes.push_back(indices[7]);
		nodes.push_back(indices[6]);
		nodes.push_back(indices[4]);
		break;
	default:
		return;
	}
}




