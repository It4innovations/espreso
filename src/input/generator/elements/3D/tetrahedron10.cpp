
#include "tetrahedron10.h"

using namespace espreso::input;

size_t Tetrahedron10::subelements = 6;
size_t Tetrahedron10::subnodes[] = { 3, 3, 3 };

void Tetrahedron10::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal size = 10;
	eslocal tetra[size];
	tetra[0] = indices[2];
	tetra[1] = indices[6];
	tetra[2] = indices[0];
	tetra[3] = indices[20];

	tetra[4] = indices[4];
	tetra[5] = indices[3];
	tetra[6] = indices[1];
	tetra[7] = indices[11];
	tetra[8] = indices[13];
	tetra[9] = indices[10];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[6];
	tetra[1] = indices[0];
	tetra[2] = indices[20];
	tetra[3] = indices[18];

	tetra[4] = indices[3];
	tetra[5] = indices[10];
	tetra[6] = indices[13];
	tetra[7] = indices[12];
	tetra[8] = indices[9];
	tetra[9] = indices[19];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[24];
	tetra[1] = indices[6];
	tetra[2] = indices[20];
	tetra[3] = indices[18];

	tetra[4] = indices[15];
	tetra[5] = indices[13];
	tetra[6] = indices[22];
	tetra[7] = indices[21];
	tetra[8] = indices[12];
	tetra[9] = indices[19];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[6];
	tetra[1] = indices[26];
	tetra[2] = indices[24];
	tetra[3] = indices[20];

	tetra[4] = indices[16];
	tetra[5] = indices[25];
	tetra[6] = indices[15];
	tetra[7] = indices[13];
	tetra[8] = indices[23];
	tetra[9] = indices[22];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[8];
	tetra[1] = indices[26];
	tetra[2] = indices[6];
	tetra[3] = indices[20];

	tetra[4] = indices[17];
	tetra[5] = indices[16];
	tetra[6] = indices[7];
	tetra[7] = indices[14];
	tetra[8] = indices[23];
	tetra[9] = indices[13];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));

	tetra[0] = indices[2];
	tetra[1] = indices[20];
	tetra[2] = indices[8];
	tetra[3] = indices[6];

	tetra[4] = indices[11];
	tetra[5] = indices[14];
	tetra[6] = indices[5];
	tetra[7] = indices[4];
	tetra[8] = indices[13];
	tetra[9] = indices[7];
	elements.push_back(new espreso::Tetrahedron10(tetra, size, params));
}

void Tetrahedron10::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal triangle1[6], triangle2[6];

	switch (face) {
	case CubeFace::X_1:
		triangle1[0] = indices[ 2];
		triangle1[1] = indices[ 8];
		triangle1[2] = indices[20];
		triangle1[3] = indices[ 5];
		triangle1[4] = indices[14];
		triangle1[5] = indices[11];

		triangle2[0] = indices[ 8];
		triangle2[1] = indices[26];
		triangle2[2] = indices[20];
		triangle2[3] = indices[17];
		triangle2[4] = indices[23];
		triangle2[5] = indices[14];
		break;
	case CubeFace::Y_1:
		triangle1[0] = indices[ 8];
		triangle1[1] = indices[ 6];
		triangle1[2] = indices[26];
		triangle1[3] = indices[ 7];
		triangle1[4] = indices[16];
		triangle1[5] = indices[17];

		triangle2[0] = indices[ 6];
		triangle2[1] = indices[24];
		triangle2[2] = indices[26];
		triangle2[3] = indices[15];
		triangle2[4] = indices[25];
		triangle2[5] = indices[16];
		break;
	case CubeFace::X_0:
		triangle1[0] = indices[ 0];
		triangle1[1] = indices[ 6];
		triangle1[2] = indices[18];
		triangle1[3] = indices[ 3];
		triangle1[4] = indices[12];
		triangle1[5] = indices[ 9];

		triangle2[0] = indices[ 6];
		triangle2[1] = indices[24];
		triangle2[2] = indices[18];
		triangle2[3] = indices[15];
		triangle2[4] = indices[21];
		triangle2[5] = indices[12];
		break;
	case CubeFace::Y_0:
		triangle1[0] = indices[ 0];
		triangle1[1] = indices[ 2];
		triangle1[2] = indices[20];
		triangle1[3] = indices[ 1];
		triangle1[4] = indices[11];
		triangle1[5] = indices[10];

		triangle2[0] = indices[ 0];
		triangle2[1] = indices[20];
		triangle2[2] = indices[18];
		triangle2[3] = indices[10];
		triangle2[4] = indices[19];
		triangle2[5] = indices[ 9];
		break;
	case CubeFace::Z_0:
		triangle1[0] = indices[ 0];
		triangle1[1] = indices[ 6];
		triangle1[2] = indices[ 2];
		triangle1[3] = indices[ 3];
		triangle1[4] = indices[ 4];
		triangle1[5] = indices[ 1];

		triangle2[0] = indices[ 6];
		triangle2[1] = indices[ 8];
		triangle2[2] = indices[ 2];
		triangle2[3] = indices[ 7];
		triangle2[4] = indices[ 5];
		triangle2[5] = indices[ 4];
		break;
	case CubeFace::Z_1:
		triangle1[0] = indices[20];
		triangle1[1] = indices[26];
		triangle1[2] = indices[24];
		triangle1[3] = indices[23];
		triangle1[4] = indices[25];
		triangle1[5] = indices[22];

		triangle2[0] = indices[20];
		triangle2[1] = indices[24];
		triangle2[2] = indices[18];
		triangle2[3] = indices[22];
		triangle2[4] = indices[21];
		triangle2[5] = indices[19];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}

	faces.push_back(new espreso::Triangle6(triangle1));
	faces.push_back(new espreso::Triangle6(triangle2));
}

void Tetrahedron10::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	switch (face) {
	case CubeFace::X_1:
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[ 5]]);
		selection.push_back(nodes[indices[ 8]]);
		selection.push_back(nodes[indices[11]]);
		selection.push_back(nodes[indices[14]]);
		selection.push_back(nodes[indices[17]]);
		selection.push_back(nodes[indices[20]]);
		selection.push_back(nodes[indices[23]]);
		selection.push_back(nodes[indices[26]]);
		break;
	case CubeFace::Y_1:
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[ 7]]);
		selection.push_back(nodes[indices[ 8]]);
		selection.push_back(nodes[indices[15]]);
		selection.push_back(nodes[indices[16]]);
		selection.push_back(nodes[indices[17]]);
		selection.push_back(nodes[indices[24]]);
		selection.push_back(nodes[indices[25]]);
		selection.push_back(nodes[indices[26]]);
		break;
	case CubeFace::X_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 3]]);
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[ 9]]);
		selection.push_back(nodes[indices[12]]);
		selection.push_back(nodes[indices[15]]);
		selection.push_back(nodes[indices[18]]);
		selection.push_back(nodes[indices[21]]);
		selection.push_back(nodes[indices[24]]);
		break;
	case CubeFace::Y_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 1]]);
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[ 9]]);
		selection.push_back(nodes[indices[10]]);
		selection.push_back(nodes[indices[11]]);
		selection.push_back(nodes[indices[18]]);
		selection.push_back(nodes[indices[19]]);
		selection.push_back(nodes[indices[20]]);
		break;
	case CubeFace::Z_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 1]]);
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[ 3]]);
		selection.push_back(nodes[indices[ 4]]);
		selection.push_back(nodes[indices[ 5]]);
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[ 7]]);
		selection.push_back(nodes[indices[ 8]]);
		break;
	case CubeFace::Z_1:
		selection.push_back(nodes[indices[18]]);
		selection.push_back(nodes[indices[19]]);
		selection.push_back(nodes[indices[20]]);
		selection.push_back(nodes[indices[21]]);
		selection.push_back(nodes[indices[22]]);
		selection.push_back(nodes[indices[23]]);
		selection.push_back(nodes[indices[24]]);
		selection.push_back(nodes[indices[25]]);
		selection.push_back(nodes[indices[26]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}

