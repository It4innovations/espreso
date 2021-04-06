
#include "hexahedron8.h"

using namespace espreso;

Hexahedron8Generator::Hexahedron8Generator()
{
	subelements = 1;
	enodes = 8;
	code = Element::CODE::HEXA8;
}

void Hexahedron8Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
	elements.push_back(indices[0]);
	elements.push_back(indices[1]);
	elements.push_back(indices[3]);
	elements.push_back(indices[2]);
	elements.push_back(indices[4]);
	elements.push_back(indices[5]);
	elements.push_back(indices[7]);
	elements.push_back(indices[6]);
}

void Hexahedron8Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
	pushNodes(elements, indices, face);
	esize.push_back(4);
	etype.push_back((int)Element::CODE::SQUARE4);
}



