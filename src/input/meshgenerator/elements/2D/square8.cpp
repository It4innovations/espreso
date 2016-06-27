
#include "square8.h"

using namespace espreso::input;

size_t Square8::subelements = 1;
size_t Square8::subnodes[] = { 1, 1, 1 };

void Square8::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal square[8];
	square[0] = indices[0];
	square[1] = indices[2];
	square[2] = indices[8];
	square[3] = indices[6];

	square[4] = indices[1];
	square[5] = indices[5];
	square[6] = indices[7];
	square[7] = indices[3];
	elements.push_back(new espreso::Square8(square, params));
}



