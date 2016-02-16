
#include "hexahedron20.h"

using namespace esinput;

size_t Hexahedron20::subelements = 1;
size_t Hexahedron20::subnodes[] = { 1, 1, 1 };

void Hexahedron20::addElements(std::vector<mesh::Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal hexa[20];
	hexa[0] = indices[2];
	hexa[1] = indices[8];
	hexa[2] = indices[6];
	hexa[3] = indices[0];
	hexa[4] = indices[20];
	hexa[5] = indices[26];
	hexa[6] = indices[24];
	hexa[7] = indices[18];

	hexa[8] = indices[5];
	hexa[9] = indices[7];
	hexa[10] = indices[3];
	hexa[11] = indices[1];
	hexa[12] = indices[23];
	hexa[13] = indices[25];
	hexa[14] = indices[21];
	hexa[15] = indices[19];
	hexa[16] = indices[11];
	hexa[17] = indices[17];
	hexa[18] = indices[15];
	hexa[19] = indices[9];
	elements.push_back(new mesh::Hexahedron20(hexa, 20, params));
}



