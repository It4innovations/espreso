
#include "utils.h"

using namespace espreso::input;

espreso::Element* AnsysUtils::createElement(eslocal *indices, eslocal n, eslocal *params)
{
	Element *e = NULL;
	if (Tetrahedron4::match(indices, n)) {
		e = new Tetrahedron4(indices, n, params);
	}
	if (Tetrahedron10::match(indices, n)) {
		e = new Tetrahedron10(indices, n, params);
	}
	if (Hexahedron8::match(indices, n)) {
		e = new Hexahedron8(indices, n, params);
	}
	if (Hexahedron20::match(indices, n)) {
		e = new Hexahedron20(indices, n, params);
	}
	if (Prisma6::match(indices, n)) {
		e = new Prisma6(indices, n, params);
	}
	if (Prisma15::match(indices, n)) {
		e = new Prisma15(indices, n, params);
	}
	if (Pyramid5::match(indices, n)) {
		e = new Pyramid5(indices, n, params);
	}
	if (Pyramid13::match(indices, n)) {
		e = new Pyramid13(indices, n, params);
	}

	if (e == NULL) {
		auto print_indices = [&] () {
			std::stringstream ss;
			for (eslocal i = 0; i < n; i++) {
				ss << indices[i] << " ";
			};
			return ss.str();
		};

		ESINFO(ERROR) << "Unknown element with indices: " << print_indices();
	}

	return e;
}


