
#include "square8.h"

#include "../../../../mesh/elements/line/line3.h"
#include "../../../../mesh/elements/plane/square8.h"

using namespace espreso::input;

size_t Square8::subelements = 1;
size_t Square8::subnodes[] = { 3, 3, 1 };

void Square8::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal square[8];
	square[0] = indices[0];
	square[1] = indices[2];
	square[2] = indices[8];
	square[3] = indices[6];

	square[4] = indices[1];
	square[5] = indices[5];
	square[6] = indices[7];
	square[7] = indices[3];
	elements.push_back(new espreso::Square8(square, params));
}

void Square8::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	ESINFO(GLOBAL_ERROR) << "Generator: plane element has no faces.";
}

void Square8::addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdge edge)
{
	eslocal line[3];
	switch (edge) {
	case CubeEdge::X_0_Z_0:
		line[0] = indices[6];
		line[1] = indices[3];
		line[2] = indices[0];
		break;
	case CubeEdge::X_1_Z_0:
		line[0] = indices[2];
		line[1] = indices[5];
		line[2] = indices[8];
		break;
	case CubeEdge::Y_0_Z_0:
		line[0] = indices[0];
		line[1] = indices[1];
		line[2] = indices[2];
		break;
	case CubeEdge::Y_1_Z_0:
		line[0] = indices[8];
		line[1] = indices[7];
		line[2] = indices[6];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown edge";
	}
	edges.push_back(new espreso::Line3(line));
}

void Square8::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeEdge edge)
{
	switch (edge) {
	case CubeEdge::X_0_Z_0:
		selection.push_back(nodes[indices[6]]);
		selection.push_back(nodes[indices[3]]);
		selection.push_back(nodes[indices[0]]);
		break;
	case CubeEdge::X_1_Z_0:
		selection.push_back(nodes[indices[2]]);
		selection.push_back(nodes[indices[5]]);
		selection.push_back(nodes[indices[8]]);
		break;
	case CubeEdge::Y_0_Z_0:
		selection.push_back(nodes[indices[0]]);
		selection.push_back(nodes[indices[1]]);
		selection.push_back(nodes[indices[2]]);
		break;
	case CubeEdge::Y_1_Z_0:
		selection.push_back(nodes[indices[8]]);
		selection.push_back(nodes[indices[7]]);
		selection.push_back(nodes[indices[6]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown edge";
	}
}

void Square8::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	ESINFO(GLOBAL_ERROR) << "Implement pickNodes for a face for SQUARE*";
}

