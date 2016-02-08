#include "prisma6.h"

using namespace mesh;

std::vector< std::vector< double> > Prisma6_rst()
{
	std::vector< std::vector<double> > rst(3, std::vector<double>(Prisma6GPCount));

	switch (Prisma6GPCount) {
	case 9: {
		double v1 = 1.0 / 6.0;
		double v2 = 4.0 / 6.0;
		double v3 = sqrt(3.0 / 5.0);
		double v4 = 0.0;
		double r[] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
		double s[] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
		double t[] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };
		rst[0].assign(r, r + 9);
		rst[1].assign(s, s + 9);
		rst[2].assign(t, t + 9);
		return rst;
	}
	default:
		ESLOG(eslog::ERROR) << "Unknown number of Prisma6 GP count.";
		exit(EXIT_FAILURE);
	}
}


std::vector< std::vector< double> > _prisma6_rst = Prisma6_rst();

std::vector<DenseMatrix> Prisma6_dN() {
	std::vector<DenseMatrix> dN(
		Prisma6GPCount,
		DenseMatrix(Point::size(), Prisma6NodesCount)
	);

	for (unsigned int i = 0; i < Prisma6GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		double r = _prisma6_rst[0][i];
		double s = _prisma6_rst[1][i];
		double t = _prisma6_rst[2][i];

		// dNr - derivation of basis function
		m(0, 0) =  t / 2.0 - 1.0 / 2.0;
		m(0, 1) = -t / 2.0 + 1.0 / 2.0;
		m(0, 2) =  0.0;
		m(0, 3) = -t / 2.0 - 1.0 / 2.0;
		m(0, 4) =  t / 2.0 + 1.0 / 2.0;
		m(0, 5) =  0;

		// dNs - derivation of basis function
		m(1, 0) =  t / 2.0 - 1.0 / 2.0;
		m(1, 1) =  0.0;
		m(1, 2) = -t / 2.0 + 1.0 / 2.0;
		m(1, 3) = -t / 2.0 - 1.0 / 2.0;
		m(1, 4) =  0.0;
		m(1, 5) =  t / 2.0 + 1.0 / 2.0;

		// dNt - derivation of basis function
		m(2, 0) =  r / 2.0 + s / 2.0 - 1.0 / 2.0;
		m(2, 1) = -r / 2.0;
		m(2, 2) =         -s / 2.0;
		m(2, 3) = -r / 2.0 - s / 2.0 + 1.0 / 2.0;
		m(2, 4) =  r / 2.0;
		m(2, 5) =          s / 2.0;
	}

	return dN;
}

std::vector<DenseMatrix> Prisma6_N() {
	std::vector<DenseMatrix> N(
		Prisma6GPCount,
		DenseMatrix(1, Prisma6NodesCount)
	);

	for (unsigned int i = 0; i < Prisma6GPCount; i++) {
		DenseMatrix &m = N[i];

		double r = _prisma6_rst[0][i];
		double s = _prisma6_rst[1][i];
		double t = _prisma6_rst[2][i];

		// basis function
		m(0, 0) = 0.5 * ((1.0 - t) * (1.0 - r - s));
		m(0, 1) = 0.5 * ((1.0 - t) * r);
		m(0, 2) = 0.5 * ((1.0 - t) * s);
		m(0, 3) = 0.5 * ((1.0 + t) * (1.0 - r - s));
		m(0, 4) = 0.5 * ((1.0 + t) * r);
		m(0, 5) = 0.5 * ((1.0 + t) * s);
	}

	return N;
}

std::vector<double> Prisma6_weight()
{
	switch (Prisma6GPCount) {
	case 9: {
		double v1 = 5.0 / 54.0;
		double v2 = 8.0 / 54.0;
		double w[9] = { v1, v1, v1, v2, v2, v2, v1, v1, v1 };
		return std::vector<double> (w, w + 9);
	}
	default:
		ESLOG(eslog::ERROR) << "Unknown number of Tatrahedron10 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Prisma6::_dN = Prisma6_dN();
std::vector<DenseMatrix> Prisma6::_N = Prisma6_N();
std::vector<double> Prisma6::_weighFactor = Prisma6_weight();

bool Prisma6::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Hexahedron8 is 3D element
	return false;
#endif

	if (n != 8) {
		return false;
	}

	if (!Element::match(indices, 2, 3)) {
		return false;
	}
	if (!Element::match(indices, 6, 7)) {
		return false;
	}

	eslocal various[6] = { 0, 1, 2, 4, 5, 6 };
	for (eslocal i = 0; i < 5; i++) {
		for (eslocal j = i + 1; j < 6; j++) {
			if (Element::match(indices, various[i], various[j])) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Prisma6::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(3);

	if (nodeIndex > 2) {
		result[0] = _indices[nodeIndex - 3];
		result[1] = _indices[(nodeIndex + 1) % 3 + 3];
		result[2] = _indices[(nodeIndex + 2) % 3 + 3];
	} else {
		result[0] = _indices[nodeIndex + 3];
		result[1] = _indices[(nodeIndex + 1) % 3];
		result[2] = _indices[(nodeIndex + 2) % 3];
	}

	return result;
}

std::vector<eslocal> Prisma6::getFace(size_t face) const
{
	// bottom
	if (face == 3) {
		std::vector<eslocal> result(3);
		result[0] = _indices[1];
		result[1] = _indices[0];
		result[2] = _indices[2];
		return result;
	}

	// top
	if (face == 4) {
		std::vector<eslocal> result(3);
		result[0] = _indices[3];
		result[1] = _indices[4];
		result[2] = _indices[5];
		return result;
	}

	//sides
	std::vector<eslocal> result(4);
	result[0] = _indices[ face              ];
	result[1] = _indices[(face + 1) % 3     ];
	result[2] = _indices[(face + 1) % 3 + 3 ];
	result[3] = _indices[ face + 3          ];
	return result;
}

Prisma6::Prisma6(const eslocal *indices, eslocal n, const eslocal *params): Element(params)
{
	switch (n) {
	case 6:
		memcpy(_indices, indices, 6 * sizeof(eslocal));
		break;
	case 8:
		_indices[0] = indices[0];
		_indices[1] = indices[1];
		_indices[2] = indices[2];
		_indices[3] = indices[4];
		_indices[4] = indices[5];
		_indices[5] = indices[6];
		break;
	default:
		ESLOG(eslog::ERROR) << "It is not possible to create Prisma6 from " << n << " elements.";
	}
}

Prisma6::Prisma6(std::ifstream &is): Element(is)
{
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * size());
}






