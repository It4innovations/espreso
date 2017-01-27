
#include "triangle3.h"

#include "../../../../mesh/elements/line/line2.h"
#include "../../../../mesh/elements/plane/triangle3.h"

using namespace espreso::input;

size_t Triangle3::subelements = 2;
size_t Triangle3::subnodes[] = { 2, 2, 1 };

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

void Triangle3::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	ESINFO(GLOBAL_ERROR) << "Generator: plane element has no faces.";
}

void Triangle3::addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdge edge)
{
	eslocal line[2];
	switch (edge) {
	case CubeEdge::X_0_Z_0:
		line[0] = indices[2];
		line[1] = indices[0];
		break;
	case CubeEdge::X_1_Z_0:
		line[0] = indices[1];
		line[1] = indices[3];
		break;
	case CubeEdge::Y_0_Z_0:
		line[0] = indices[0];
		line[1] = indices[1];
		break;
	case CubeEdge::Y_1_Z_0:
		line[0] = indices[3];
		line[1] = indices[2];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown edge";
	}
	edges.push_back(new espreso::Line2(line));
}

void Triangle3::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeEdge edge)
{
	switch (edge) {
	case CubeEdge::X_0_Z_0:
		selection.push_back(nodes[indices[2]]);
		selection.push_back(nodes[indices[0]]);
		break;
	case CubeEdge::X_1_Z_0:
		selection.push_back(nodes[indices[1]]);
		selection.push_back(nodes[indices[3]]);
		break;
	case CubeEdge::Y_0_Z_0:
		selection.push_back(nodes[indices[0]]);
		selection.push_back(nodes[indices[1]]);
		break;
	case CubeEdge::Y_1_Z_0:
		selection.push_back(nodes[indices[3]]);
		selection.push_back(nodes[indices[2]]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown edge";
	}
}

void Triangle3::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	ESINFO(GLOBAL_ERROR) << "Implement pickNodes for a face for TRIANGLE3";
}

