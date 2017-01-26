
#include "triangle6.h"

#include "../../../../mesh/elements/line/line3.h"
#include "../../../../mesh/elements/plane/triangle6.h"

using namespace espreso::input;

size_t Triangle6::subelements = 2;
size_t Triangle6::subnodes[] = { 3, 3, 1 };

void Triangle6::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal triangle[6];

	triangle[0] = indices[0];
	triangle[1] = indices[2];
	triangle[2] = indices[8];

	triangle[3] = indices[1];
	triangle[4] = indices[5];
	triangle[5] = indices[4];
	elements.push_back(new espreso::Triangle6(triangle, params));

	triangle[0] = indices[0];
	triangle[1] = indices[8];
	triangle[2] = indices[6];

	triangle[3] = indices[4];
	triangle[4] = indices[7];
	triangle[5] = indices[3];
	elements.push_back(new espreso::Triangle6(triangle, params));
}

void Triangle6::addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdge edge)
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

void Triangle6::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeEdge edge)
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

