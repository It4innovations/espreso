
#include "triangle3.h"

using namespace espreso::input;

size_t Triangle3::subelements = 2;
size_t Triangle3::subnodes[] = { 0, 0, 0 };

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



