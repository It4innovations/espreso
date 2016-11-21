
#include "hexahedron8.h"

using namespace espreso::input;

size_t Hexahedron8::subelements = 1;
size_t Hexahedron8::subnodes[] = { 0, 0, 0 };

void Hexahedron8::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal hexa[8];
	hexa[0] = indices[0];
	hexa[1] = indices[1];
	hexa[2] = indices[3];
	hexa[3] = indices[2];
	hexa[4] = indices[4];
	hexa[5] = indices[5];
	hexa[6] = indices[7];
	hexa[7] = indices[6];
	elements.push_back(new espreso::Hexahedron8(hexa, 8, params));
}

void Hexahedron8::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal square[4];

	switch (face) {
	case CubeFace::X_1:
		square[0] = indices[1];
		square[1] = indices[3];
		square[2] = indices[7];
		square[3] = indices[5];
		break;
	case CubeFace::Y_1:
		square[0] = indices[3];
		square[1] = indices[2];
		square[2] = indices[6];
		square[3] = indices[7];
		break;
	case CubeFace::X_0:
		square[0] = indices[2];
		square[1] = indices[0];
		square[2] = indices[4];
		square[3] = indices[6];
		break;
	case CubeFace::Y_0:
		square[0] = indices[0];
		square[1] = indices[1];
		square[2] = indices[5];
		square[3] = indices[4];
		break;
	case CubeFace::Z_0:
		square[0] = indices[0];
		square[1] = indices[2];
		square[2] = indices[3];
		square[3] = indices[1];
		break;
	case CubeFace::Z_1:
		square[0] = indices[5];
		square[1] = indices[7];
		square[2] = indices[6];
		square[3] = indices[4];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}

	faces.push_back(new espreso::Square4(square));
}

void Hexahedron8::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	switch (face) {
	case CubeFace::X_1:
		selection.push_back(nodes[indices[1]]);
		selection.push_back(nodes[indices[3]]);
		selection.push_back(nodes[indices[7]]);
		selection.push_back(nodes[indices[5]]);
		break;
	case CubeFace::Y_1:
		selection.push_back(nodes[indices[3]]);
		selection.push_back(nodes[indices[2]]);
		selection.push_back(nodes[indices[6]]);
		selection.push_back(nodes[indices[7]]);
		break;
	case CubeFace::X_0:
		selection.push_back(nodes[indices[2]]);
		selection.push_back(nodes[indices[0]]);
		selection.push_back(nodes[indices[4]]);
		selection.push_back(nodes[indices[6]]);
		break;
	case CubeFace::Y_0:
		selection.push_back(nodes[indices[0]]);
		selection.push_back(nodes[indices[1]]);
		selection.push_back(nodes[indices[5]]);
		selection.push_back(nodes[indices[4]]);
		break;
	case CubeFace::Z_0:
		selection.push_back(nodes[indices[0]]);
		selection.push_back(nodes[indices[2]]);
		selection.push_back(nodes[indices[3]]);
		selection.push_back(nodes[indices[1]]);
		break;
	case CubeFace::Z_1:
		selection.push_back(nodes[indices[5]]);
		selection.push_back(nodes[indices[7]]);
		selection.push_back(nodes[indices[6]]);
		selection.push_back(nodes[indices[4]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}



