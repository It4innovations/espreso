
#include "utils.h"

using namespace espreso::input;

static std::string print_indices(eslocal *indices, eslocal n)
{
	std::stringstream ss;
	for (eslocal i = 0; i < n; i++) {
		ss << indices[i] << " ";
	};
	return ss.str();
};

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
		ESINFO(ERROR) << "Unknown element with indices: " << print_indices(indices, n);
	}

	return e;
}

espreso::Element* AnsysUtils::createElement(eslocal *indices, eslocal n, eslocal *params, int eType)
{
	switch (eType) {
	case 185:
		if (Hexahedron8::match(indices, n)) {
			return new Hexahedron8(indices, n, params);
		}
		if (Tetrahedron4::match(indices, n)) {
			return new Tetrahedron4(indices, n, params);
		}
		if (Prisma6::match(indices, n)) {
			return new Prisma6(indices, n, params);
		}
		if (Pyramid5::match(indices, n)) {
			return new Pyramid5(indices, n, params);
		}
		ESINFO(ERROR) << "Unknown element with indices: " << print_indices(indices, n);
	case 90:
	case 186:
		if (Hexahedron20::match(indices, n)) {
			return new Hexahedron20(indices, n, params);
		}
		if (Tetrahedron10::match(indices, n)) {
			return new Tetrahedron10(indices, n, params);
		}
		if (Prisma15::match(indices, n)) {
			return new Prisma15(indices, n, params);
		}
		if (Pyramid13::match(indices, n)) {
			return new Pyramid13(indices, n, params);
		}
		ESINFO(ERROR) << "Unknown element with indices: " << print_indices(indices, n);
	case 187:
		return new Tetrahedron10(indices, n, params);
	case 152:
		if (Square4::match(indices, n)) {
			return new Square4(indices, params);
		}
		if (Triangle3::match(indices, n)) {
			return new Triangle3(indices, params);
		}
	case 154:
		if (Square8::match(indices, n)) {
			return new Square8(indices, params);
		}
		if (Triangle6::match(indices, n)) {
			return new Triangle6(indices, params);
		}
		ESINFO(ERROR) << "Unknown element with indices: " << print_indices(indices, n);
	default:
		ESINFO(ERROR) << "Unknown element type: " << eType;
	}

	return NULL;
}


