
#include "pyramid13.h"

using namespace espreso::input;

size_t Pyramid13::subelements = 6;
size_t Pyramid13::subnodes[] = { 3, 3, 3 };

void Pyramid13::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal pyramid[13];
	pyramid[0] = indices[100];
	pyramid[1] = indices[104];
	pyramid[2] = indices[4];
	pyramid[3] = indices[0];
	pyramid[4] = indices[62];

	pyramid[5] = indices[102];
	pyramid[6] = indices[54];
	pyramid[7] = indices[2];
	pyramid[8] = indices[50];
	pyramid[9] = indices[81];
	pyramid[10] = indices[83];
	pyramid[11] = indices[33];
	pyramid[12] = indices[31];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[104];
	pyramid[1] = indices[124];
	pyramid[2] = indices[24];
	pyramid[3] = indices[4];

	pyramid[5] = indices[114];
	pyramid[6] = indices[74];
	pyramid[7] = indices[14];
	pyramid[8] = indices[54];
	pyramid[9] = indices[83];
	pyramid[10] = indices[93];
	pyramid[11] = indices[43];
	pyramid[12] = indices[33];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[124];
	pyramid[1] = indices[120];
	pyramid[2] = indices[20];
	pyramid[3] = indices[24];

	pyramid[5] = indices[122];
	pyramid[6] = indices[70];
	pyramid[7] = indices[22];
	pyramid[8] = indices[74];
	pyramid[9] = indices[93];
	pyramid[10] = indices[91];
	pyramid[11] = indices[41];
	pyramid[12] = indices[43];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[120];
	pyramid[1] = indices[100];
	pyramid[2] = indices[0];
	pyramid[3] = indices[20];

	pyramid[5] = indices[110];
	pyramid[6] = indices[50];
	pyramid[7] = indices[10];
	pyramid[8] = indices[70];
	pyramid[9] = indices[91];
	pyramid[10] = indices[81];
	pyramid[11] = indices[31];
	pyramid[12] = indices[41];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[100];
	pyramid[1] = indices[120];
	pyramid[2] = indices[124];
	pyramid[3] = indices[104];

	pyramid[5] = indices[110];
	pyramid[6] = indices[122];
	pyramid[7] = indices[114];
	pyramid[8] = indices[102];
	pyramid[9] = indices[81];
	pyramid[10] = indices[91];
	pyramid[11] = indices[93];
	pyramid[12] = indices[83];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[4];
	pyramid[1] = indices[24];
	pyramid[2] = indices[20];
	pyramid[3] = indices[0];

	pyramid[5] = indices[14];
	pyramid[6] = indices[22];
	pyramid[7] = indices[10];
	pyramid[8] = indices[2];
	pyramid[9] = indices[33];
	pyramid[10] = indices[43];
	pyramid[11] = indices[41];
	pyramid[12] = indices[31];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));
}

void Pyramid13::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal square[8];

	switch (face) {
	case CubeFace::X_1:
		square[0] = indices[  4];
		square[1] = indices[ 24];
		square[2] = indices[124];
		square[3] = indices[104];
		square[4] = indices[ 14];
		square[5] = indices[ 74];
		square[6] = indices[114];
		square[7] = indices[ 54];
		break;
	case CubeFace::Y_1:
		square[0] = indices[ 24];
		square[1] = indices[ 20];
		square[2] = indices[120];
		square[3] = indices[124];
		square[4] = indices[ 22];
		square[5] = indices[ 70];
		square[6] = indices[122];
		square[7] = indices[ 74];
		break;
	case CubeFace::X_0:
		square[0] = indices[ 20];
		square[1] = indices[  0];
		square[2] = indices[100];
		square[3] = indices[120];
		square[4] = indices[ 10];
		square[5] = indices[ 50];
		square[6] = indices[110];
		square[7] = indices[ 70];
		break;
	case CubeFace::Y_0:
		square[0] = indices[  0];
		square[1] = indices[  4];
		square[2] = indices[104];
		square[3] = indices[100];
		square[4] = indices[  2];
		square[5] = indices[ 54];
		square[6] = indices[102];
		square[7] = indices[ 50];
		break;
	case CubeFace::Z_0:
		square[0] = indices[  0];
		square[1] = indices[ 20];
		square[2] = indices[ 24];
		square[3] = indices[  4];
		square[4] = indices[ 10];
		square[5] = indices[ 22];
		square[6] = indices[ 14];
		square[7] = indices[  2];
		break;
	case CubeFace::Z_1:
		square[0] = indices[104];
		square[1] = indices[124];
		square[2] = indices[120];
		square[3] = indices[100];
		square[4] = indices[114];
		square[5] = indices[122];
		square[6] = indices[110];
		square[7] = indices[102];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}

	faces.push_back(new espreso::Square8(square));
}

void Pyramid13::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	switch (face) {
	case CubeFace::X_1:
		selection.push_back(nodes[indices[  4]]);
		selection.push_back(nodes[indices[ 14]]);
		selection.push_back(nodes[indices[ 24]]);
		selection.push_back(nodes[indices[ 54]]);
		selection.push_back(nodes[indices[ 74]]);
		selection.push_back(nodes[indices[104]]);
		selection.push_back(nodes[indices[114]]);
		selection.push_back(nodes[indices[124]]);
		break;
	case CubeFace::Y_1:
		selection.push_back(nodes[indices[ 20]]);
		selection.push_back(nodes[indices[ 22]]);
		selection.push_back(nodes[indices[ 24]]);
		selection.push_back(nodes[indices[ 70]]);
		selection.push_back(nodes[indices[ 74]]);
		selection.push_back(nodes[indices[120]]);
		selection.push_back(nodes[indices[122]]);
		selection.push_back(nodes[indices[124]]);
		break;
	case CubeFace::X_0:
		selection.push_back(nodes[indices[  0]]);
		selection.push_back(nodes[indices[ 10]]);
		selection.push_back(nodes[indices[ 20]]);
		selection.push_back(nodes[indices[ 50]]);
		selection.push_back(nodes[indices[ 70]]);
		selection.push_back(nodes[indices[100]]);
		selection.push_back(nodes[indices[110]]);
		selection.push_back(nodes[indices[120]]);
		break;
	case CubeFace::Y_0:
		selection.push_back(nodes[indices[  0]]);
		selection.push_back(nodes[indices[  2]]);
		selection.push_back(nodes[indices[  4]]);
		selection.push_back(nodes[indices[ 50]]);
		selection.push_back(nodes[indices[ 54]]);
		selection.push_back(nodes[indices[100]]);
		selection.push_back(nodes[indices[102]]);
		selection.push_back(nodes[indices[104]]);
		break;
	case CubeFace::Z_0:
		selection.push_back(nodes[indices[  0]]);
		selection.push_back(nodes[indices[  2]]);
		selection.push_back(nodes[indices[  4]]);
		selection.push_back(nodes[indices[ 10]]);
		selection.push_back(nodes[indices[ 14]]);
		selection.push_back(nodes[indices[ 20]]);
		selection.push_back(nodes[indices[ 22]]);
		selection.push_back(nodes[indices[ 24]]);
		break;
	case CubeFace::Z_1:
		selection.push_back(nodes[indices[100]]);
		selection.push_back(nodes[indices[102]]);
		selection.push_back(nodes[indices[104]]);
		selection.push_back(nodes[indices[110]]);
		selection.push_back(nodes[indices[114]]);
		selection.push_back(nodes[indices[120]]);
		selection.push_back(nodes[indices[122]]);
		selection.push_back(nodes[indices[124]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}

