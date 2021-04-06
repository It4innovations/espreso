
#include "square4.h"

using namespace espreso;

Square4Generator::Square4Generator()
{
	subelements = 1;
	enodes = 4;
	code = Element::CODE::SQUARE4;
}

void Square4Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
	elements.push_back(indices[0]);
	elements.push_back(indices[1]);
	elements.push_back(indices[3]);
	elements.push_back(indices[2]);
}



