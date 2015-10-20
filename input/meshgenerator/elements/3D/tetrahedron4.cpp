
#include "tetrahedron4.h"

using namespace esinput;

size_t Tetrahedron4::subelements = 6;
size_t Tetrahedron4::subnodes[] = { 0, 0, 0 };

void Tetrahedron4::addElements(std::vector<mesh::Element*> &elements, const eslocal indices[])
{
	eslocal tetra[5];
	tetra[0] = indices[0];
	tetra[1] = indices[3];
	tetra[2] = indices[2];
	tetra[4] = indices[4];
	elements.push_back(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[3];
	tetra[1] = indices[2];
	tetra[2] = indices[4];
	tetra[4] = indices[6];
	elements.push_back(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[7];
	tetra[1] = indices[3];
	tetra[2] = indices[4];
	tetra[4] = indices[6];
	elements.push_back(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[3];
	tetra[1] = indices[5];
	tetra[2] = indices[7];
	tetra[4] = indices[4];
	elements.push_back(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[1];
	tetra[1] = indices[5];
	tetra[2] = indices[3];
	tetra[4] = indices[4];
	elements.push_back(new mesh::Tetrahedron4(tetra));

	tetra[0] = indices[0];
	tetra[1] = indices[4];
	tetra[2] = indices[1];
	tetra[4] = indices[3];
	elements.push_back(new mesh::Tetrahedron4(tetra));
}



