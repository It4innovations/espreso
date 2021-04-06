
#include "prisma15.h"

using namespace espreso;

Prisma15Generator::Prisma15Generator()
{
	subelements = 2;
	enodes = 15;
	code = Element::CODE::PRISMA15;
}

void Prisma15Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
	elements.push_back(indices[ 0]);
	elements.push_back(indices[ 2]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[18]);
	elements.push_back(indices[20]);
	elements.push_back(indices[24]);

	elements.push_back(indices[ 1]);
	elements.push_back(indices[ 4]);
	elements.push_back(indices[ 3]);
	elements.push_back(indices[19]);
	elements.push_back(indices[22]);
	elements.push_back(indices[21]);
	elements.push_back(indices[ 9]);
	elements.push_back(indices[11]);
	elements.push_back(indices[15]);


	elements.push_back(indices[ 2]);
	elements.push_back(indices[ 8]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[20]);
	elements.push_back(indices[26]);
	elements.push_back(indices[24]);

	elements.push_back(indices[ 5]);
	elements.push_back(indices[ 7]);
	elements.push_back(indices[ 4]);
	elements.push_back(indices[23]);
	elements.push_back(indices[25]);
	elements.push_back(indices[22]);
	elements.push_back(indices[11]);
	elements.push_back(indices[17]);
	elements.push_back(indices[15]);
}

void Prisma15Generator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
	switch (face) {
	case CubeFace::X_0:
	case CubeFace::X_1:
	case CubeFace::Y_0:
	case CubeFace::Y_1:
		pushSquareNodes(nodes, indices, face);
		break;
	case CubeFace::Z_0:
	case CubeFace::Z_1:
		pushTriangleNodes(nodes, indices, face);
		break;
	default:
		break;
	}
}

void Prisma15Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
	pushNodes(elements, indices, face);
	switch (face) {
	case CubeFace::X_0:
	case CubeFace::X_1:
	case CubeFace::Y_0:
	case CubeFace::Y_1:
		esize.push_back(8);
		etype.push_back((int)Element::CODE::SQUARE8);
		break;
	case CubeFace::Z_0:
	case CubeFace::Z_1:
		esize.push_back(6);
		esize.push_back(6);
		etype.push_back((int)Element::CODE::TRIANGLE6);
		etype.push_back((int)Element::CODE::TRIANGLE6);
		break;
	default:
		break;
	}
}



