#include "hexahedron8.h"

using namespace mesh;

std::vector<DenseMatrix> Hexa_dN() {
	std::vector<DenseMatrix> dN(
		Hexahedron8GPCount,
		DenseMatrix(Point::size(), Hexahedron8NodesCount)
	);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Hexahedron8GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNr - derivation of basis function
		m(0, 0) = 0.125 * (-(1 - s) * (1 - t));
		m(0, 1) = 0.125 * ( (1 - s) * (1 - t));
		m(0, 2) = 0.125 * ( (1 + s) * (1 - t));
		m(0, 3) = 0.125 * (-(1 + s) * (1 - t));
		m(0, 4) = 0.125 * (-(1 - s) * (1 + t));
		m(0, 5) = 0.125 * ( (1 - s) * (1 + t));
		m(0, 6) = 0.125 * ( (1 + s) * (1 + t));
		m(0, 7) = 0.125 * (-(1 + s) * (1 + t));

		// dNs - derivation of basis function
		m(1, 0)  = 0.125 * (-(1 - r) * (1 - t));
		m(1, 1)  = 0.125 * (-(1 + r) * (1 - t));
		m(1, 2) = 0.125 * ( (1 + r) * (1 - t));
		m(1, 3) = 0.125 * ( (1 - r) * (1 - t));
		m(1, 4) = 0.125 * (-(1 - r) * (1 + t));
		m(1, 5) = 0.125 * (-(1 + r) * (1 + t));
		m(1, 6) = 0.125 * ( (1 + r) * (1 + t));
		m(1, 7) = 0.125 * ( (1 - r) * (1 + t));

		// dNt - derivation of basis function
		m(2, 0) = 0.125 * (-(1 - r) * (1 - s));
		m(2, 1) = 0.125 * (-(1 + r) * (1 - s));
		m(2, 2) = 0.125 * (-(1 + r) * (1 + s));
		m(2, 3) = 0.125 * (-(1 - r) * (1 + s));
		m(2, 4) = 0.125 * ( (1 - r) * (1 - s));
		m(2, 5) = 0.125 * ( (1 + r) * (1 - s));
		m(2, 6) = 0.125 * ( (1 + r) * (1 + s));
		m(2, 7) = 0.125 * ( (1 - r) * (1 + s));
	}

	return dN;
}

std::vector<DenseMatrix> Hexa_N() {
	std::vector<DenseMatrix> N(
		Hexahedron8GPCount,
		DenseMatrix(1, Hexahedron8NodesCount)
	);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Hexahedron8GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		N[i](0, 0) = 0.125 * (1 - r) * (1 - s) * (1 - t);
		N[i](0, 1) = 0.125 * (r + 1) * (1 - s) * (1 - t);
		N[i](0, 2) = 0.125 * (r + 1) * (s + 1) * (1 - t);
		N[i](0, 3) = 0.125 * (1 - r) * (s + 1) * (1 - t);
		N[i](0, 4) = 0.125 * (1 - r) * (1 - s) * (t + 1);
		N[i](0, 5) = 0.125 * (r + 1) * (1 - s) * (t + 1);
		N[i](0, 6) = 0.125 * (r + 1) * (s + 1) * (t + 1);
		N[i](0, 7) = 0.125 * (1 - r) * (s + 1) * (t + 1);
	}

	return N;
}

std::vector<DenseMatrix> Hexahedron8::_dN = Hexa_dN();
std::vector<DenseMatrix> Hexahedron8::_N = Hexa_N();
std::vector<double> Hexahedron8::_weighFactor(Hexahedron8NodesCount, 1);

bool Hexahedron8::match(esint *indices, esint n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Hexahedron8 is 3D element
	return false;
#endif

	if (n != 8) {
		return false;
	}

	for (esint i = 0; i < 7; i++) {
		for (esint j = i + 1; j < 8; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<esint> Hexahedron8::getNeighbours(size_t nodeIndex) const
{
	std::vector<esint> result(3);

	if (nodeIndex > 3) {
		result[0] = _indices[nodeIndex - 4];
		result[1] = _indices[(nodeIndex + 1) % 4 + 4];
		result[2] = _indices[(nodeIndex + 3) % 4 + 4];
	} else {
		result[0] = _indices[nodeIndex + 4];
		result[1] = _indices[(nodeIndex + 1) % 4];
		result[2] = _indices[(nodeIndex + 3) % 4];
	}

	return result;
}

std::vector<esint> Hexahedron8::getFace(size_t face) const
{
	std::vector<esint> result(4);

	// bottom
	if (face == 4) {
		result[0] = _indices[0];
		result[1] = _indices[3];
		result[2] = _indices[2];
		result[3] = _indices[1];
		return result;
	}

	// top
	if (face == 5) {
		result[0] = _indices[4];
		result[1] = _indices[5];
		result[2] = _indices[6];
		result[3] = _indices[7];
		return result;
	}

	//sides
	result[0] = _indices[ face              ];
	result[1] = _indices[(face + 1) % 4     ];
	result[2] = _indices[(face + 1) % 4 + 4 ];
	result[3] = _indices[ face + 4          ];
	return result;
}

Hexahedron8::Hexahedron8(esint *indices)
{
	memcpy(_indices, indices, Hexahedron8NodesCount * sizeof(esint));
}


