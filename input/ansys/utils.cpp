
#include "utils.h"

using namespace esinput;

mesh::Element* AnsysUtils::createElement(eslocal *indices, eslocal n, eslocal *params)
{
	mesh::Element *e = NULL;
	if (mesh::Tetrahedron4::match(indices, n)) {
		e = new mesh::Tetrahedron4(indices, n, params);
	}
	if (mesh::Tetrahedron10::match(indices, n)) {
		e = new mesh::Tetrahedron10(indices, n, params);
	}
	if (mesh::Hexahedron8::match(indices, n)) {
		e = new mesh::Hexahedron8(indices, n, params);
	}
	if (mesh::Hexahedron20::match(indices, n)) {
		e = new mesh::Hexahedron20(indices, n, params);
	}
	if (mesh::Prisma6::match(indices, n)) {
		e = new mesh::Prisma6(indices, n, params);
	}
	if (mesh::Prisma15::match(indices, n)) {
		e = new mesh::Prisma15(indices, n, params);
	}
	if (mesh::Pyramid5::match(indices, n)) {
		e = new mesh::Pyramid5(indices, n, params);
	}
	if (mesh::Pyramid13::match(indices, n)) {
		e = new mesh::Pyramid13(indices, n, params);
	}

	if (e == NULL) {
		auto print_indices = [&] () {
			std::stringstream ss;
			for (eslocal i = 0; i < n; i++) {
				ss << indices[i] << " ";
			};
			return ss.str();
		};

		ESINFO(eslog::ERROR) << "Unknown element with indices: " << print_indices();
	}

	return e;
}


