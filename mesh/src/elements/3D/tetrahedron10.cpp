#include "tetrahedron10.h"

std::vector<std::vector<double> > Tetra_dN()
{
	std::vector<std::vector<double> > dN(Tetrahedron10GPCount);

	int dN_length  = 40;

	double CsQ_scale = 2.236067977499790;

	double diag = (5.0 + 3.0 * CsQ_scale) / 20.0;
	double rest = (5.0 - CsQ_scale) / 20.0;

	for (unsigned int i = 0; i < Tetrahedron10GPCount; i++) {
		double r = (i == 3) ? diag : rest;
		double s = (i == 2) ? diag : rest;
		double t = (i == 1) ? diag : rest;
		double q = (i == 0) ? diag : rest;

		dN[i].assign(dN_length, 0);

		// dNq = [ 4 * q - 1, 0, 0, 0, 4 * r, 0, 4 * s, 4 * t, 0, 0 ];
		dN[i][0] = 4.0 * q - 1;
		dN[i][4] = 4.0 * r;
		dN[i][6] = 4.0 * s;
		dN[i][7] = 4.0 * t;

		// dNr = [0, 4 * r - 1, 0, 0, 4 * q, 4 * s, 0, 0, 4 * t, 0 ];
		dN[i][11] = 4.0 * r - 1.0;
		dN[i][14] = 4.0 * q;
		dN[i][15] = 4.0 * s;
		dN[i][18] = 4.0 * t;

		// dNs = [ 0, 0, 4 * s - 1, 0, 0, 4 * r,  4 * q, 0, 0, 4 * t ];
		dN[i][22] = 4.0 * s - 1.0;
		dN[i][25] = 4.0 * r;
		dN[i][26] = 4.0 * q;
		dN[i][29] = 4.0 * t;

		// dNt = [ 0, 0, 0, 4 * t - 1, 0, 0, 0, 4 * q, 4 * r, 4 * s];
		dN[i][33] = 4.0 * t - 1.0;
		dN[i][37] = 4.0 * q;
		dN[i][38] = 4.0 * r;
		dN[i][39] = 4.0 * s;
	}

	return dN;
}

std::vector<std::vector<double> > Tetra_N() {
	std::vector<std::vector<double> > N(Tetrahedron10GPCount);

	double CsQ_scale = 0.577350269189626;

	double diag = (5.0 + 3.0 * CsQ_scale) / 20.0;
	double rest = (5.0 - CsQ_scale) / 20.0;

	for (unsigned int i = 0; i < Tetrahedron10GPCount; i++) {
		double r = (i == 3) ? diag : rest;
		double s = (i == 2) ? diag : rest;
		double t = (i == 1) ? diag : rest;
		double q = (i == 0) ? diag : rest;

		N[i].resize(10);
		N[i][0] = (2.0 * q - 1.0) * q;
		N[i][1] = (2.0 * r - 1.0) * r;
		N[i][2] = (2.0 * s - 1.0) * s;
		N[i][3] = (2.0 * t - 1.0) * t;
		N[i][4] = 4.0 * q * r;
		N[i][5] = 4.0 * r * s;
		N[i][6] = 4.0 * q * s;
		N[i][7] = 4.0 * q * t;
		N[i][8] = 4.0 * r * t;
		N[i][9] = 4.0 * s * t;
	}

	return N;
}

std::vector<std::vector<double> > Tetrahedron10::_dN = Tetra_dN();
std::vector<std::vector<double> > Tetrahedron10::_N = Tetra_N();
std::vector<double> Tetrahedron10::_weighFactor(Tetrahedron10NodesCount, 1.0 / 24.0);

bool Tetrahedron10::match(idx_t *indices, idx_t n) {

#ifndef D3
	// Tetrahedron10 is 3D element
	return false;
#endif

	if (n != 20) {
		return false;
	}

	if (!Element::match(indices, 2, 3)) {
		return false;
	}
	if (!Element::match(indices, 2, 10)) {
		return false;
	}
	if (!Element::match(indices, 19, 20)) {
		return false;
	}
	if (!Element::match(indices, 4, 5)) {
		return false;
	}
	if (!Element::match(indices, 4, 6)) {
		return false;
	}
	if (!Element::match(indices, 4, 7)) {
		return false;
	}
	if (!Element::match(indices, 4, 12)) {
		return false;
	}
	if (!Element::match(indices, 4, 13)) {
		return false;
	}
	if (!Element::match(indices, 4, 14)) {
		return false;
	}
	if (!Element::match(indices, 4, 15)) {
		return false;
	}

	idx_t various[10] = { 0, 1, 2, 4, 8, 9 ,11, 16, 17, 18 };
	for (idx_t i = 0; i < 9; i++) {
		for (idx_t j = i + 1; j < 10; j++) {
			if (Element::match(indices, various[i], various[j])) {
				return false;
			}
		}
	}

	return true;
}

std::vector<idx_t> Tetrahedron10::getNeighbours(size_t nodeIndex) const
{
	std::vector<idx_t> result;
	if (nodeIndex > 3) {
		result.resize(2);
	} else {
		result.resize(3);
	}

	switch (nodeIndex) {
	case 0: {
		result[0] = _indices[4];
		result[1] = _indices[6];
		result[2] = _indices[7];
		return result;
	}
	case 1: {
		result[0] = _indices[4];
		result[1] = _indices[5];
		result[2] = _indices[8];
		return result;
	}
	case 2: {
		result[0] = _indices[5];
		result[1] = _indices[6];
		result[2] = _indices[9];
		return result;
	}
	case 3: {
		result[0] = _indices[7];
		result[1] = _indices[8];
		result[2] = _indices[9];
		return result;
	}
	case 4: {
		result[0] = _indices[0];
		result[1] = _indices[1];
		return result;
	}
	case 5: {
		result[0] = _indices[1];
		result[1] = _indices[2];
		return result;
	}
	case 6: {
		result[0] = _indices[0];
		result[1] = _indices[2];
		return result;
	}
	case 7: {
		result[0] = _indices[0];
		result[1] = _indices[3];
		return result;
	}
	case 8: {
		result[0] = _indices[1];
		result[1] = _indices[3];
		return result;
	}
	case 9: {
		result[0] = _indices[2];
		result[1] = _indices[3];
		return result;
	}
	}
	return result;
}

std::vector<idx_t> Tetrahedron10::getFace(size_t face) const
{
	std::vector<idx_t> result(3);
	result[0] = (face < 3) ? _indices[0] : _indices[1];
	result[1] = (face < 2) ? _indices[1] : _indices[2];
	result[2] = (face < 1) ? _indices[2] : _indices[3];
	return result;
}

Tetrahedron10::Tetrahedron10(idx_t *indices)
{
	_indices[0] = indices[0];
	_indices[1] = indices[1];
	_indices[2] = indices[2];
	_indices[3] = indices[4];
	_indices[4] = indices[8];
	_indices[5] = indices[9];
	_indices[6] = indices[11];
	_indices[7] = indices[16];
	_indices[8] = indices[17];
	_indices[9] = indices[18];
}



