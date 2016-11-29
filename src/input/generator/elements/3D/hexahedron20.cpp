
#include "hexahedron20.h"

using namespace espreso::input;

size_t Hexahedron20::subelements = 1;
size_t Hexahedron20::subnodes[] = { 3, 3, 3 };

void Hexahedron20::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal hexa[20];
	hexa[0] = indices[2];
	hexa[1] = indices[8];
	hexa[2] = indices[6];
	hexa[3] = indices[0];
	hexa[4] = indices[20];
	hexa[5] = indices[26];
	hexa[6] = indices[24];
	hexa[7] = indices[18];

	hexa[8] = indices[5];
	hexa[9] = indices[7];
	hexa[10] = indices[3];
	hexa[11] = indices[1];
	hexa[12] = indices[23];
	hexa[13] = indices[25];
	hexa[14] = indices[21];
	hexa[15] = indices[19];
	hexa[16] = indices[11];
	hexa[17] = indices[17];
	hexa[18] = indices[15];
	hexa[19] = indices[9];
	elements.push_back(new espreso::Hexahedron20(hexa, 20, params));
}

void Hexahedron20::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal square[8];

	switch (face) {
	case CubeFace::X_1:
		square[0] = indices[ 2];
		square[1] = indices[ 8];
		square[2] = indices[26];
		square[3] = indices[20];

		square[4] = indices[ 5];
		square[5] = indices[17];
		square[6] = indices[23];
		square[7] = indices[11];
		break;
	case CubeFace::Y_1:
		square[0] = indices[ 8];
		square[1] = indices[ 6];
		square[2] = indices[24];
		square[3] = indices[26];

		square[4] = indices[ 7];
		square[5] = indices[15];
		square[6] = indices[25];
		square[7] = indices[17];
		break;
	case CubeFace::X_0:
		square[0] = indices[ 6];
		square[1] = indices[ 0];
		square[2] = indices[18];
		square[3] = indices[24];

		square[4] = indices[ 3];
		square[5] = indices[ 9];
		square[6] = indices[21];
		square[7] = indices[15];
		break;
	case CubeFace::Y_0:
		square[0] = indices[ 0];
		square[1] = indices[ 2];
		square[2] = indices[20];
		square[3] = indices[18];

		square[4] = indices[ 1];
		square[5] = indices[11];
		square[6] = indices[19];
		square[7] = indices[ 9];
		break;
	case CubeFace::Z_0:
		square[0] = indices[ 0];
		square[1] = indices[ 6];
		square[2] = indices[ 8];
		square[3] = indices[ 2];

		square[4] = indices[ 3];
		square[5] = indices[ 7];
		square[6] = indices[ 5];
		square[7] = indices[ 1];
		break;
	case CubeFace::Z_1:
		square[0] = indices[20];
		square[1] = indices[26];
		square[2] = indices[24];
		square[3] = indices[18];

		square[4] = indices[23];
		square[5] = indices[25];
		square[6] = indices[21];
		square[7] = indices[19];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}

	faces.push_back(new espreso::Square8(square));
}

void Hexahedron20::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	switch (face) {
	case CubeFace::X_1:
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[ 8]]);
		selection.push_back(nodes[indices[26]]);
		selection.push_back(nodes[indices[20]]);

		selection.push_back(nodes[indices[ 5]]);
		selection.push_back(nodes[indices[17]]);
		selection.push_back(nodes[indices[23]]);
		selection.push_back(nodes[indices[11]]);
		break;
	case CubeFace::Y_1:
		selection.push_back(nodes[indices[ 8]]);
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[24]]);
		selection.push_back(nodes[indices[26]]);

		selection.push_back(nodes[indices[ 7]]);
		selection.push_back(nodes[indices[15]]);
		selection.push_back(nodes[indices[25]]);
		selection.push_back(nodes[indices[17]]);
		break;
	case CubeFace::X_0:
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[18]]);
		selection.push_back(nodes[indices[24]]);

		selection.push_back(nodes[indices[ 3]]);
		selection.push_back(nodes[indices[ 9]]);
		selection.push_back(nodes[indices[21]]);
		selection.push_back(nodes[indices[15]]);
		break;
	case CubeFace::Y_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[20]]);
		selection.push_back(nodes[indices[18]]);

		selection.push_back(nodes[indices[ 1]]);
		selection.push_back(nodes[indices[11]]);
		selection.push_back(nodes[indices[19]]);
		selection.push_back(nodes[indices[ 9]]);
		break;
	case CubeFace::Z_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[ 8]]);
		selection.push_back(nodes[indices[ 2]]);

		selection.push_back(nodes[indices[ 3]]);
		selection.push_back(nodes[indices[ 7]]);
		selection.push_back(nodes[indices[ 5]]);
		selection.push_back(nodes[indices[ 1]]);
		break;
	case CubeFace::Z_1:
		selection.push_back(nodes[indices[20]]);
		selection.push_back(nodes[indices[26]]);
		selection.push_back(nodes[indices[24]]);
		selection.push_back(nodes[indices[18]]);

		selection.push_back(nodes[indices[23]]);
		selection.push_back(nodes[indices[25]]);
		selection.push_back(nodes[indices[21]]);
		selection.push_back(nodes[indices[19]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}


