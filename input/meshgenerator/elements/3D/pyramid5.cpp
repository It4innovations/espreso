
#include "pyramid5.h"

using namespace esinput;

size_t Pyramid5::subelements = 6;
size_t Pyramid5::subnodes[] = { 1, 1, 1 };

void Pyramid5::addElements(std::vector<mesh::Element*> &elements, const eslocal indices[])
{
	eslocal pyramid[8];
	pyramid[0] = indices[18];
	pyramid[1] = indices[20];
	pyramid[2] = indices[2];
	pyramid[3] = indices[0];
	pyramid[4] = indices[13];
	pyramid[5] = indices[13];
	pyramid[6] = indices[13];
	pyramid[7] = indices[13];
	elements.push_back(new mesh::Pyramid5(pyramid));

	pyramid[0] = indices[20];
	pyramid[1] = indices[26];
	pyramid[2] = indices[8];
	pyramid[3] = indices[2];
	elements.push_back(new mesh::Pyramid5(pyramid));

	pyramid[0] = indices[26];
	pyramid[1] = indices[24];
	pyramid[2] = indices[6];
	pyramid[3] = indices[8];
	elements.push_back(new mesh::Pyramid5(pyramid));

	pyramid[0] = indices[24];
	pyramid[1] = indices[18];
	pyramid[2] = indices[0];
	pyramid[3] = indices[6];
	elements.push_back(new mesh::Pyramid5(pyramid));

	pyramid[0] = indices[18];
	pyramid[1] = indices[24];
	pyramid[2] = indices[26];
	pyramid[3] = indices[20];
	elements.push_back(new mesh::Pyramid5(pyramid));

	pyramid[0] = indices[2];
	pyramid[1] = indices[8];
	pyramid[2] = indices[6];
	pyramid[3] = indices[0];
	elements.push_back(new mesh::Pyramid5(pyramid));
}



