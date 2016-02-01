
#include "prisma15.h"

using namespace esinput;

size_t Prisma15::subelements = 2;
size_t Prisma15::subnodes[] = { 1, 1, 1 };

void Prisma15::addElements(std::vector<mesh::Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal prisma[20];
	prisma[0] = indices[0];
	prisma[1] = indices[2];
	prisma[2] = indices[8];
	prisma[3] = indices[8];
	prisma[4] = indices[18];
	prisma[5] = indices[20];
	prisma[6] = indices[26];
	prisma[7] = indices[26];

	prisma[8] = indices[1];
	prisma[9] = indices[5];
	prisma[10] = indices[8];
	prisma[11] = indices[4];
	prisma[12] = indices[19];
	prisma[13] = indices[23];
	prisma[14] = indices[26];
	prisma[15] = indices[22];
	prisma[16] = indices[9];
	prisma[17] = indices[11];
	prisma[18] = indices[17];
	prisma[19] = indices[17];
	elements.push_back(new mesh::Prisma15(prisma, params));

	prisma[0] = indices[0];
	prisma[1] = indices[8];
	prisma[2] = indices[6];
	prisma[3] = indices[6];
	prisma[4] = indices[18];
	prisma[5] = indices[26];
	prisma[6] = indices[24];
	prisma[7] = indices[24];

	prisma[8] = indices[4];
	prisma[9] = indices[7];
	prisma[10] = indices[6];
	prisma[11] = indices[3];
	prisma[12] = indices[22];
	prisma[13] = indices[25];
	prisma[14] = indices[24];
	prisma[15] = indices[21];
	prisma[16] = indices[9];
	prisma[17] = indices[17];
	prisma[18] = indices[15];
	prisma[19] = indices[15];
	elements.push_back(new mesh::Prisma15(prisma, params));
}



