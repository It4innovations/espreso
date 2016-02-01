
#include "hexahedron8.h"

using namespace esinput;

size_t Hexahedron8::subelements = 1;
size_t Hexahedron8::subnodes[] = { 0, 0, 0 };

void Hexahedron8::addElements(std::vector<mesh::Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal hexa[8];
	hexa[0] = indices[0];
	hexa[1] = indices[1];
	hexa[2] = indices[3];
	hexa[3] = indices[2];
	hexa[4] = indices[4];
	hexa[5] = indices[5];
	hexa[6] = indices[7];
	hexa[7] = indices[6];
	elements.push_back(new mesh::Hexahedron8(hexa, params));
}



