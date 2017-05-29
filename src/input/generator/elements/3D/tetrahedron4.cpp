
#include "tetrahedron4.h"
#include "hexahedron8.h"

#include "../../../../mesh/elements/plane/triangle3.h"
#include "../../../../mesh/elements/volume/tetrahedron4.h"

using namespace espreso::input;

size_t Tetrahedron4::subelements = 6;
size_t Tetrahedron4::subnodes[] = { 2, 2, 2 };

void Tetrahedron4::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal tetra[4];
	tetra[0] = indices[0];
	tetra[1] = indices[3];
	tetra[2] = indices[2];
	tetra[3] = indices[4];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[3];
	tetra[1] = indices[2];
	tetra[2] = indices[4];
	tetra[3] = indices[6];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[7];
	tetra[1] = indices[3];
	tetra[2] = indices[4];
	tetra[3] = indices[6];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[3];
	tetra[1] = indices[5];
	tetra[2] = indices[7];
	tetra[3] = indices[4];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[1];
	tetra[1] = indices[5];
	tetra[2] = indices[3];
	tetra[3] = indices[4];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));

	tetra[0] = indices[0];
	tetra[1] = indices[4];
	tetra[2] = indices[1];
	tetra[3] = indices[3];
	elements.push_back(new espreso::Tetrahedron4(tetra, 4, params));
}

void Tetrahedron4::addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdge edge)
{
	ESINFO(GLOBAL_ERROR) << "Implement addEdges for TETRAHEDRON4";
}

void Tetrahedron4::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal triangle1[3], triangle2[3];

	switch (face) {
	case CubeFace::X_1:
		triangle1[0] = indices[1];
		triangle1[1] = indices[3];
		triangle1[2] = indices[5];

		triangle2[0] = indices[3];
		triangle2[1] = indices[7];
		triangle2[2] = indices[5];
		break;
	case CubeFace::Y_1:
		triangle1[0] = indices[3];
		triangle1[1] = indices[2];
		triangle1[2] = indices[6];

		triangle2[0] = indices[3];
		triangle2[1] = indices[6];
		triangle2[2] = indices[7];
		break;
	case CubeFace::X_0:
		triangle1[0] = indices[2];
		triangle1[1] = indices[0];
		triangle1[2] = indices[4];

		triangle2[0] = indices[2];
		triangle2[1] = indices[4];
		triangle2[2] = indices[6];
		break;
	case CubeFace::Y_0:
		triangle1[0] = indices[0];
		triangle1[1] = indices[1];
		triangle1[2] = indices[4];

		triangle2[0] = indices[1];
		triangle2[1] = indices[5];
		triangle2[2] = indices[4];
		break;
	case CubeFace::Z_0:
		triangle1[0] = indices[0];
		triangle1[1] = indices[2];
		triangle1[2] = indices[3];

		triangle2[0] = indices[0];
		triangle2[1] = indices[3];
		triangle2[2] = indices[1];
		break;
	case CubeFace::Z_1:
		triangle1[0] = indices[5];
		triangle1[1] = indices[7];
		triangle1[2] = indices[4];

		triangle2[0] = indices[7];
		triangle2[1] = indices[6];
		triangle2[2] = indices[4];
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}

	faces.push_back(new espreso::Triangle3(triangle1));
	faces.push_back(new espreso::Triangle3(triangle2));
}

void Tetrahedron4::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeEdge edge)
{
	ESINFO(GLOBAL_ERROR) << "Implement pickNodes for an edge for HEXA20";
}

void Tetrahedron4::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	Hexahedron8::pickNodes(nodes, selection, indices, face);
}
