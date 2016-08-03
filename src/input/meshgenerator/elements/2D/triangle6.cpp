
#include "triangle6.h"

using namespace espreso::input;

size_t Triangle6::subelements = 2;
size_t Triangle6::subnodes[] = { 1, 1, 1 };

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



