
#include "pyramid5.h"

#include "../../../../mesh/elements/plane/square4.h"
#include "../../../../mesh/elements/volume/pyramid5.h"

using namespace espreso::input;

size_t Pyramid5::subelements = 6;
size_t Pyramid5::subnodes[] = { 3, 3, 3 };

void Pyramid5::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal pyramid[5];
	pyramid[0] = indices[18];
	pyramid[1] = indices[20];
	pyramid[2] = indices[2];
	pyramid[3] = indices[0];
	pyramid[4] = indices[13];
	elements.push_back(new espreso::Pyramid5(pyramid, 5, params));

	pyramid[0] = indices[20];
	pyramid[1] = indices[26];
	pyramid[2] = indices[8];
	pyramid[3] = indices[2];
	elements.push_back(new espreso::Pyramid5(pyramid, 5, params));

	pyramid[0] = indices[26];
	pyramid[1] = indices[24];
	pyramid[2] = indices[6];
	pyramid[3] = indices[8];
	elements.push_back(new espreso::Pyramid5(pyramid, 5, params));

	pyramid[0] = indices[24];
	pyramid[1] = indices[18];
	pyramid[2] = indices[0];
	pyramid[3] = indices[6];
	elements.push_back(new espreso::Pyramid5(pyramid, 5, params));

	pyramid[0] = indices[18];
	pyramid[1] = indices[24];
	pyramid[2] = indices[26];
	pyramid[3] = indices[20];
	elements.push_back(new espreso::Pyramid5(pyramid, 5, params));

	pyramid[0] = indices[2];
	pyramid[1] = indices[8];
	pyramid[2] = indices[6];
	pyramid[3] = indices[0];
	elements.push_back(new espreso::Pyramid5(pyramid, 5, params));
}


void Pyramid5::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal square[4];

	switch (face) {
	case CubeFace::X_1:
		square[0] = indices[ 2];
		square[1] = indices[ 8];
		square[2] = indices[26];
		square[3] = indices[20];
		break;
	case CubeFace::Y_1:
		square[0] = indices[ 8];
		square[1] = indices[ 6];
		square[2] = indices[24];
		square[3] = indices[26];
		break;
	case CubeFace::X_0:
		square[0] = indices[ 6];
		square[1] = indices[ 0];
		square[2] = indices[18];
		square[3] = indices[24];
		break;
	case CubeFace::Y_0:
		square[0] = indices[ 0];
		square[1] = indices[ 2];
		square[2] = indices[20];
		square[3] = indices[18];
		break;
	case CubeFace::Z_0:
		square[0] = indices[ 0];
		square[1] = indices[ 6];
		square[2] = indices[ 8];
		square[3] = indices[ 2];
		break;
	case CubeFace::Z_1:
		square[0] = indices[20];
		square[1] = indices[26];
		square[2] = indices[24];
		square[3] = indices[18];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}

	faces.push_back(new espreso::Square4(square));
}

void Pyramid5::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	switch (face) {
	case CubeFace::X_1:
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[ 8]]);
		selection.push_back(nodes[indices[20]]);
		selection.push_back(nodes[indices[26]]);
		break;
	case CubeFace::Y_1:
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[ 8]]);
		selection.push_back(nodes[indices[24]]);
		selection.push_back(nodes[indices[26]]);
		break;
	case CubeFace::X_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[18]]);
		selection.push_back(nodes[indices[24]]);
		break;
	case CubeFace::Y_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[18]]);
		selection.push_back(nodes[indices[20]]);
		break;
	case CubeFace::Z_0:
		selection.push_back(nodes[indices[ 0]]);
		selection.push_back(nodes[indices[ 2]]);
		selection.push_back(nodes[indices[ 6]]);
		selection.push_back(nodes[indices[ 8]]);
		break;
	case CubeFace::Z_1:
		selection.push_back(nodes[indices[18]]);
		selection.push_back(nodes[indices[20]]);
		selection.push_back(nodes[indices[24]]);
		selection.push_back(nodes[indices[26]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}
