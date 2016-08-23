
#include "square4.h"

using namespace espreso::input;

size_t Square4::subelements = 1;
size_t Square4::subnodes[] = { 0, 0, 0 };

void Square4::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal square[4];
	square[0] = indices[0];
	square[1] = indices[1];
	square[2] = indices[3];
	square[3] = indices[2];
	elements.push_back(new espreso::Square4(square, params));
}

void Square4::addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdges edge)
{
	eslocal line[2];
	switch (edge) {
	case CubeEdges::X_0_Z_0:
		line[0] = indices[2];
		line[1] = indices[0];
		break;
	case CubeEdges::X_1_Z_0:
		line[0] = indices[1];
		line[1] = indices[3];
		break;
	case CubeEdges::Y_0_Z_0:
		line[0] = indices[0];
		line[1] = indices[1];
		break;
	case CubeEdges::Y_1_Z_0:
		line[0] = indices[3];
		line[1] = indices[2];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown edge";
	}
	edges.push_back(new espreso::Line2(line));
}



