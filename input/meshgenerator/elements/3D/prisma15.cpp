
#include "prisma15.h"

using namespace espreso::input;

size_t Prisma15::subelements = 2;
size_t Prisma15::subnodes[] = { 1, 1, 1 };

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



