
#include "square4.h"

using namespace espreso::input;

size_t Square4::subelements = 1;
size_t Square4::subnodes[] = { 0, 0, 0 };

void Square4::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal square[4];
	square[0] = indices[0];
	square[1] = indices[1];
	square[2] = indices[3];
	square[3] = indices[2];
	elements.push_back(new espreso::Square4(square, params));
}



