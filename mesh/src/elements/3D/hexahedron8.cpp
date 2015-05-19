#include "hexahedron8.h"

std::vector<std::vector<double> > Hexa_dN() {
	std::vector<std::vector<double> > dN(Hexahedron8GPCount);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Hexahedron8GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		///dN contains [dNr, dNs, dNt]
		dN[i].resize(Point::size() * Hexahedron8NodesCount);

		// dNr - derivation of basis function
		dN[i][0] = 0.125 * (-(1 - s) * (1 - t));
		dN[i][1] = 0.125 * ( (1 - s) * (1 - t));
		dN[i][2] = 0.125 * ( (1 + s) * (1 - t));
		dN[i][3] = 0.125 * (-(1 + s) * (1 - t));
		dN[i][4] = 0.125 * (-(1 - s) * (1 + t));
		dN[i][5] = 0.125 * ( (1 - s) * (1 + t));
		dN[i][6] = 0.125 * ( (1 + s) * (1 + t));
		dN[i][7] = 0.125 * (-(1 + s) * (1 + t));

		// dNs - derivation of basis function
		dN[i][8]  = 0.125 * (-(1 - r) * (1 - t));
		dN[i][9]  = 0.125 * (-(1 + r) * (1 - t));
		dN[i][10] = 0.125 * ( (1 + r) * (1 - t));
		dN[i][11] = 0.125 * ( (1 - r) * (1 - t));
		dN[i][12] = 0.125 * (-(1 - r) * (1 + t));
		dN[i][13] = 0.125 * (-(1 + r) * (1 + t));
		dN[i][14] = 0.125 * ( (1 + r) * (1 + t));
		dN[i][15] = 0.125 * ( (1 - r) * (1 + t));

		// dNt - derivation of basis function
		dN[i][16] = 0.125 * (-(1 - r) * (1 - s));
		dN[i][17] = 0.125 * (-(1 + r) * (1 - s));
		dN[i][18] = 0.125 * (-(1 + r) * (1 + s));
		dN[i][19] = 0.125 * (-(1 - r) * (1 + s));
		dN[i][20] = 0.125 * ( (1 - r) * (1 - s));
		dN[i][21] = 0.125 * ( (1 + r) * (1 - s));
		dN[i][22] = 0.125 * ( (1 + r) * (1 + s));
		dN[i][23] = 0.125 * ( (1 - r) * (1 + s));
	}

	return dN;
}

std::vector<std::vector<double> > Hexa_N() {
	std::vector<std::vector<double> > N(Hexahedron8GPCount);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Hexahedron8GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		N[i].resize(Hexahedron8NodesCount);
		N[i][0] = 0.125 * (1 - r) * (1 - s) * (1 - t);
		N[i][1] = 0.125 * (r + 1) * (1 - s) * (1 - t);
		N[i][2] = 0.125 * (r + 1) * (s + 1) * (1 - t);
		N[i][3] = 0.125 * (1 - r) * (s + 1) * (1 - t);
		N[i][4] = 0.125 * (1 - r) * (1 - s) * (t + 1);
		N[i][5] = 0.125 * (r + 1) * (1 - s) * (t + 1);
		N[i][6] = 0.125 * (r + 1) * (s + 1) * (t + 1);
		N[i][7] = 0.125 * (1 - r) * (s + 1) * (t + 1);
	}

	return N;
}

std::vector<std::vector<double> > Hexahedron8::_dN = Hexa_dN();
std::vector<std::vector<double> > Hexahedron8::_N = Hexa_N();
std::vector<double> Hexahedron8::_weighFactor(Hexahedron8NodesCount, 1);

bool Hexahedron8::match(idx_t *indices, idx_t n) {

#ifndef D3
	// Hexahedron8 is 3D element
	return false;
#endif

	if (n != 8) {
		return false;
	}

	for (idx_t i = 0; i < 7; i++) {
		for (idx_t j = i + 1; j < 8; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<idx_t> Hexahedron8::getNeighbours(size_t nodeIndex) const
{
	std::vector<idx_t> result(3);

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

std::vector<idx_t> Hexahedron8::getFace(size_t face) const
{
	std::vector<idx_t> result(4);

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

Hexahedron8::Hexahedron8(idx_t *indices)
{
	memcpy(_indices, indices, Hexahedron8NodesCount * sizeof(idx_t));
}


