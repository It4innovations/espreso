
#include "prisma6.h"
#include "hexahedron8.h"

using namespace espreso::input;

size_t Prisma6::subelements = 2;
size_t Prisma6::subnodes[] = { 2, 2, 2 };

void Prisma6::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal prisma[6];
	prisma[0] = indices[0];
	prisma[1] = indices[1];
	prisma[2] = indices[3];
	prisma[3] = indices[4];
	prisma[4] = indices[5];
	prisma[5] = indices[7];
	elements.push_back(new espreso::Prisma6(prisma, 6, params));

	prisma[0] = indices[0];
	prisma[1] = indices[3];
	prisma[2] = indices[2];
	prisma[3] = indices[4];
	prisma[4] = indices[7];
	prisma[5] = indices[6];
	elements.push_back(new espreso::Prisma6(prisma, 6, params));
}


void Prisma6::addFaces(std::vector<Element*> &faces, const eslocal indices[], CubeFace face)
{
	eslocal triangle1[3], triangle2[3];

	switch (face) {
	case CubeFace::X_1:
	case CubeFace::Y_1:
	case CubeFace::X_0:
	case CubeFace::Y_0:
		Hexahedron8::addFaces(faces, indices, face);
		break;
	case CubeFace::Z_0:
		triangle1[0] = indices[0];
		triangle1[1] = indices[2];
		triangle1[2] = indices[3];
		faces.push_back(new espreso::Triangle3(triangle1));

		triangle2[0] = indices[0];
		triangle2[1] = indices[3];
		triangle2[2] = indices[1];
		faces.push_back(new espreso::Triangle3(triangle2));
		break;
	case CubeFace::Z_1:
		triangle1[0] = indices[4];
		triangle1[1] = indices[5];
		triangle1[2] = indices[7];
		faces.push_back(new espreso::Triangle3(triangle1));

		triangle2[0] = indices[4];
		triangle2[1] = indices[7];
		triangle2[2] = indices[6];
		faces.push_back(new espreso::Triangle3(triangle2));
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Incorrect face";
	}
}

void Prisma6::pickNodes(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const eslocal indices[], CubeFace face)
{
	Hexahedron8::pickNodes(nodes, selection, indices, face);
}
