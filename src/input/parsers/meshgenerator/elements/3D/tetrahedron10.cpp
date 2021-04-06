
#include "tetrahedron10.h"

using namespace espreso;

Tetrahedron10Generator::Tetrahedron10Generator()
{
	subelements = 6;
	enodes = 10;
	code = Element::CODE::TETRA10;
}

void Tetrahedron10Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
	elements.push_back(indices[ 2]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[ 0]);
	elements.push_back(indices[20]);

	elements.push_back(indices[ 4]);
	elements.push_back(indices[ 3]);
	elements.push_back(indices[ 1]);
	elements.push_back(indices[11]);
	elements.push_back(indices[13]);
	elements.push_back(indices[10]);


	elements.push_back(indices[ 6]);
	elements.push_back(indices[ 0]);
	elements.push_back(indices[20]);
	elements.push_back(indices[18]);

	elements.push_back(indices[ 3]);
	elements.push_back(indices[10]);
	elements.push_back(indices[13]);
	elements.push_back(indices[12]);
	elements.push_back(indices[ 9]);
	elements.push_back(indices[19]);


	elements.push_back(indices[24]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[20]);
	elements.push_back(indices[18]);

	elements.push_back(indices[15]);
	elements.push_back(indices[13]);
	elements.push_back(indices[22]);
	elements.push_back(indices[21]);
	elements.push_back(indices[12]);
	elements.push_back(indices[19]);


	elements.push_back(indices[ 6]);
	elements.push_back(indices[26]);
	elements.push_back(indices[24]);
	elements.push_back(indices[20]);

	elements.push_back(indices[16]);
	elements.push_back(indices[25]);
	elements.push_back(indices[15]);
	elements.push_back(indices[13]);
	elements.push_back(indices[23]);
	elements.push_back(indices[22]);


	elements.push_back(indices[ 8]);
	elements.push_back(indices[26]);
	elements.push_back(indices[ 6]);
	elements.push_back(indices[20]);

	elements.push_back(indices[17]);
	elements.push_back(indices[16]);
	elements.push_back(indices[ 7]);
	elements.push_back(indices[14]);
	elements.push_back(indices[23]);
	elements.push_back(indices[13]);


	elements.push_back(indices[ 2]);
	elements.push_back(indices[20]);
	elements.push_back(indices[ 8]);
	elements.push_back(indices[ 6]);

	elements.push_back(indices[11]);
	elements.push_back(indices[14]);
	elements.push_back(indices[ 5]);
	elements.push_back(indices[ 4]);
	elements.push_back(indices[13]);
	elements.push_back(indices[ 7]);
}

void Tetrahedron10Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
	pushTriangleNodes(elements, indices, face);
	esize.push_back(6);
	esize.push_back(6);
	etype.push_back((int)Element::CODE::TRIANGLE6);
	etype.push_back((int)Element::CODE::TRIANGLE6);
}



