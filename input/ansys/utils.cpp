
#include "utils.h"

using namespace esinput;

mesh::Element* AnsysUtils::createElement(eslocal *indices, eslocal n, eslocal *params)
{
	mesh::Element *e = NULL;
	if (mesh::Tetrahedron4::match(indices, n)) {
		e = new mesh::Tetrahedron4(indices, params);
	}
	if (mesh::Tetrahedron10::match(indices, n)) {
		e = new mesh::Tetrahedron10(indices, n, params);
	}
	if (mesh::Hexahedron8::match(indices, n)) {
		e = new mesh::Hexahedron8(indices, params);
	}
	if (mesh::Hexahedron20::match(indices, n)) {
		e = new mesh::Hexahedron20(indices, params);
	}
	if (mesh::Prisma6::match(indices, n)) {
		e = new mesh::Prisma6(indices, params);
	}
	if (mesh::Prisma15::match(indices, n)) {
		e = new mesh::Prisma15(indices, params);
	}
	if (mesh::Pyramid5::match(indices, n)) {
		e = new mesh::Pyramid5(indices, params);
	}
	if (mesh::Pyramid13::match(indices, n)) {
		e = new mesh::Pyramid13(indices, params);
	}

	if (e == NULL) {
		std::cerr << "Unknown element with indices: ";
		for (eslocal i = 0; i < n; i++) {
			std::cerr << indices[i] << " ";
		}
		std::cerr << "\n";
		exit(EXIT_FAILURE);
	}

	return e;
}


