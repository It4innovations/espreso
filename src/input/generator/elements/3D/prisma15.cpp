
#include "prisma15.h"
#include "hexahedron20.h"
#include "tetrahedron10.h"

#include "../../../../mesh/elements/plane/triangle6.h"
#include "../../../../mesh/elements/volume/prisma15.h"

using namespace espreso::input;

size_t Prisma15::subelements = 2;
size_t Prisma15::subnodes[] = { 3, 3, 3 };

void Prisma15::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal prisma[15];
	prisma[0] = indices[0];
	prisma[1] = indices[2];
	prisma[2] = indices[8];
	prisma[3] = indices[18];
	prisma[4] = indices[20];
	prisma[5] = indices[26];

	prisma[6] = indices[1];
	prisma[7] = indices[5];
	prisma[8] = indices[4];
	prisma[9] = indices[19];
	prisma[10] = indices[23];
	prisma[11] = indices[22];
	prisma[12] = indices[9];
	prisma[13] = indices[11];
	prisma[14] = indices[17];
	elements.push_back(new espreso::Prisma15(prisma, 15, params));

	prisma[0] = indices[0];
	prisma[1] = indices[8];
	prisma[2] = indices[6];
	prisma[3] = indices[18];
	prisma[4] = indices[26];
	prisma[5] = indices[24];;

	prisma[6] = indices[4];
	prisma[7] = indices[7];
	prisma[8] = indices[3];
	prisma[9] = indices[22];
	prisma[10] = indices[25];
	prisma[11] = indices[21];
	prisma[12] = indices[9];
	prisma[13] = indices[17];
	prisma[14] = indices[15];
	elements.push_back(new espreso::Prisma15(prisma, 15, params));
}

void Prisma15::addEdges(std::vector<Element*> &edges, const eslocal indices[], CubeEdge edge)
{
	ESINFO(GLOBAL_ERROR) << "Implement addEdges for PRISMA15";
}

void Prisma15::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal triangle1[6], triangle2[6];

	switch (face) {
	case CubeFace::X_1:
	case CubeFace::Y_1:
	case CubeFace::X_0:
	case CubeFace::Y_0:
		Hexahedron20::addFaces(faces, indices, face);
		break;
	case CubeFace::Z_1:
		triangle1[0] = indices[20];
		triangle1[1] = indices[26];
		triangle1[2] = indices[18];
		triangle1[3] = indices[23];
		triangle1[4] = indices[22];
		triangle1[5] = indices[19];
		faces.push_back(new Triangle6(triangle1));

		triangle2[0] = indices[26];
		triangle2[1] = indices[24];
		triangle2[2] = indices[18];
		triangle2[3] = indices[25];
		triangle2[4] = indices[21];
		triangle2[5] = indices[22];
		faces.push_back(new Triangle6(triangle2));
		break;
	case CubeFace::Z_0:
		triangle1[0] = indices[ 0];
		triangle1[1] = indices[ 6];
		triangle1[2] = indices[ 8];
		triangle1[3] = indices[ 3];
		triangle1[4] = indices[ 7];
		triangle1[5] = indices[ 4];
		faces.push_back(new Triangle6(triangle1));

		triangle2[0] = indices[ 0];
		triangle2[1] = indices[ 8];
		triangle2[2] = indices[ 2];
		triangle2[3] = indices[ 4];
		triangle2[4] = indices[ 5];
		triangle2[5] = indices[ 1];
		faces.push_back(new Triangle6(triangle2));
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}

void Prisma15::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeEdge edge)
{
	ESINFO(GLOBAL_ERROR) << "Implement pickNodes for an edge for HEXA20";
}

void Prisma15::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	switch (face) {
	case CubeFace::X_1:
	case CubeFace::Y_1:
	case CubeFace::X_0:
	case CubeFace::Y_0:
		Hexahedron20::pickNodes(nodes, selection, indices, face);
		break;
	case CubeFace::Z_0:
	case CubeFace::Z_1:
		Tetrahedron10::pickNodes(nodes, selection, indices, face);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}
