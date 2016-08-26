
#include "triangle3.h"

using namespace espreso::input;

size_t Triangle3::subelements = 2;
size_t Triangle3::subnodes[] = { 0, 0, 0 };

void Triangle3::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal triangle[3];

	triangle[0] = indices[0];
	triangle[1] = indices[1];
	triangle[2] = indices[3];
	elements.push_back(new espreso::Triangle3(triangle, params));

	triangle[0] = indices[0];
	triangle[1] = indices[3];
	triangle[2] = indices[2];
	elements.push_back(new espreso::Triangle3(triangle, params));
}

void Triangle3::addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdges edge)
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

void Triangle3::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeEdges edge)
{
	switch (edge) {
	case CubeEdges::X_0_Z_0:
		selection.push_back(nodes[indices[2]]);
		selection.push_back(nodes[indices[0]]);
		break;
	case CubeEdges::X_1_Z_0:
		selection.push_back(nodes[indices[1]]);
		selection.push_back(nodes[indices[3]]);
		break;
	case CubeEdges::Y_0_Z_0:
		selection.push_back(nodes[indices[0]]);
		selection.push_back(nodes[indices[1]]);
		break;
	case CubeEdges::Y_1_Z_0:
		selection.push_back(nodes[indices[3]]);
		selection.push_back(nodes[indices[2]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown edge";
	}
}

