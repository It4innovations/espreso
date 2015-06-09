
#include "hexahedron20.h"

using namespace mesh;

std::vector< std::vector< double> > rst()
{
	std::vector< std::vector<double> > rst(3, std::vector<double>(Hexahedron20GPCount));

	switch (Hexahedron20GPCount) {
	case 8: {
		double v = 0.577350269189625953;
		double r[] = {  v,  v,  v,  v, -v, -v, -v, -v };
		double s[] = { -v, -v,  v,  v, -v, -v,  v,  v };
		double t[] = { -v,  v, -v,  v, -v,  v, -v,  v };
		rst[0].assign(r, r + 8);
		rst[1].assign(s, s + 8);
		rst[2].assign(t, t + 8);
		return rst;
	}
	case 14: {
		double v1 = 0.758786910639329015;
		double v2 = 0.795822425754222018;
		double v3 = 0;
		double r[] = {
			-v1,  v1,  v1, -v1, -v1,  v1,  v1, -v1,
			 v3,  v3,  v2,  v3, -v2,  v3
		};
		double s[] = {
			-v1, -v1,  v1,  v1, -v1, -v1,  v1,  v1,
			 v3, -v2,  v3,  v2,  v3,  v3
		};
		double t[] = {
			-v1, -v1, -v1, -v1,  v1,  v1,  v1,  v1,
			-v2,  v3,  v3,  v3,  v3,  v2
		};
		rst[0].assign(r, r + 14);
		rst[1].assign(s, s + 14);
		rst[2].assign(t, t + 14);
		return rst;
	}
	default:
		std::cerr << "Unknown number of Hexahedron20 GP count\n";
		exit(EXIT_FAILURE);
	}
}


std::vector< std::vector< double> > _rst = rst();

std::vector<DenseMatrix> Hexa20_dN()
{
	std::vector<DenseMatrix> dN(
		Hexahedron20GPCount,
		DenseMatrix(Point::size(), Hexahedron20NodesCount)
	);

	for (unsigned int i = 0; i < Hexahedron20GPCount; i++) {
		DenseMatrix &m = dN[i];

		double r = _rst[0][i];
		double s = _rst[1][i];
		double t = _rst[2][i];


		// dNr - derivation of basis function
		m(0, 0) =  ((s - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(0, 1) =  ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0 - ((s - 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0;
		m(0, 2) = -((s + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(0, 3) = -((s + 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(0, 4) = -((s - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(0, 5) = -((s - 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(0, 6) =  ((s + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(0, 7) =  ((s + 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 + ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(0, 8)  = -(r * (s - 1.0) * (t - 1.0)) / 2.0;
		m(0, 9)  =  ((pow(s, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(0, 10) =  (r * (s + 1.0) * (t - 1.0)) / 2.0;
		m(0, 11) = -((pow(s, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(0, 12) =  (r * (s - 1.0) * (t + 1.0)) / 2.0;
		m(0, 13) = -((pow(s, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(0, 14) = -(r * (s + 1.0) * (t + 1.0)) / 2.0;
		m(0, 15) =  ((pow(s, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(0, 16) = -((pow(t, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(0, 17) =  ((pow(t, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(0, 18) = -((pow(t, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(0, 19) =  ((pow(t, 2.0) - 1.0) * (s + 1.0)) / 4.0;


		// dNs - derivation of basis function
		m(1, 0) =  ((r - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(1, 1) = -((r + 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(1, 2) = -((r + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(1, 3) =  ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r - 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0;
		m(1, 4) = -((r - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(1, 5) =  ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r + 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0;
		m(1, 6) =  ((r + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(1, 7) =  ((r - 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(1, 8)  = -((pow(r, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(1, 9)  =  (s * (r + 1.0) * (t - 1.0)) / 2.0;
		m(1, 10) =  ((pow(r, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(1, 11) = -(s * (r - 1.0) * (t - 1.0)) / 2.0;
		m(1, 12) =  ((pow(r, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(1, 13) = -(s * (r + 1.0) * (t + 1.0)) / 2.0;
		m(1, 14) = -((pow(r, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(1, 15) =  (s * (r - 1.0) * (t + 1.0)) / 2.0;
		m(1, 16) = -((pow(t, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(1, 17) =  ((pow(t, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(1, 18) = -((pow(t, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(1, 19) =  ((pow(t, 2.0) - 1.0) * (r - 1.0)) / 4.0;

		// dNt - derivation of basis function
		m(2, 0) =  ((r - 1.0) * (s - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(2, 1) = -((r + 1.0) * (s - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(2, 2) =  ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r + 1.0) * (s + 1.0) * (r + s - t - 2.0)) / 8.0;
		m(2, 3) = -((r - 1.0) * (s + 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(2, 4) =  ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r - 1.0) * (s - 1.0) * (r + s - t + 2.0)) / 8.0;
		m(2, 5) = -((r + 1.0) * (s - 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(2, 6) =  ((r + 1.0) * (s + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(2, 7) =  ((r - 1.0) * (s + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(2, 8)  = -((pow(r, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(2, 9)  =  ((pow(s, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(2, 10) =  ((pow(r, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(2, 11) = -((pow(s, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(2, 12) =  ((pow(r, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(2, 13) = -((pow(s, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(2, 14) = -((pow(r, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(2, 15) =  ((pow(s, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(2, 16) = -(t * (r - 1.0) * (s - 1.0)) / 2.0;
		m(2, 17) =  (t * (r + 1.0) * (s - 1.0)) / 2.0;
		m(2, 18) = -(t * (r + 1.0) * (s + 1.0)) / 2.0;
		m(2, 19) =  (t * (r - 1.0) * (s + 1.0)) / 2.0;
	}

	return dN;
}




std::vector<DenseMatrix> Hexa20_N() {
	std::vector<DenseMatrix> N(
		Hexahedron20GPCount,
		DenseMatrix(1, Hexahedron20NodesCount)
	);

	for (unsigned int i = 0; i < Hexahedron20GPCount; i++) {
		DenseMatrix &m = N[i];

		double r = _rst[0][i];
		double s = _rst[1][i];
		double t = _rst[2][i];

		// basis function
		m(0, 0) = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 - t) * (-r - s - t - 2.0));
		m(0, 1) = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 - t) * ( r - s - t - 2.0));
		m(0, 2) = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 - t) * ( r + s - t - 2.0));
		m(0, 3) = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 - t) * (-r + s - t - 2.0));
		m(0, 4) = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 + t) * (-r - s + t - 2.0));
		m(0, 5) = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 + t) * ( r - s + t - 2.0));
		m(0, 6) = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 + t) * ( r + s + t - 2.0));
		m(0, 7) = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 + t) * (-r + s + t - 2.0));

		m(0, 8) =  0.25 * ((1.0 - pow(r, 2.0)) * (1.0 - s) * (1.0 - t));
		m(0, 9) =  0.25 * ((1.0 + r) * (1.0 - pow(s, 2.0)) * (1.0 - t));
		m(0, 10) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 + s) * (1.0 - t));
		m(0, 11) = 0.25 * ((1.0 - r) * (1.0 - pow(s, 2.0)) * (1.0 - t));
		m(0, 12) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 - s) * (1.0 + t));
		m(0, 13) = 0.25 * ((1.0 + r) * (1.0 - pow(s, 2.0)) * (1.0 + t));
		m(0, 14) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 + s) * (1.0 + t));
		m(0, 15) = 0.25 * ((1.0 - r) * (1.0 - pow(s, 2.0)) * (1.0 + t));
		m(0, 16) = 0.25 * ((1.0 - r) * (1.0 - s) * (1.0 - pow(t, 2.0)));
		m(0, 17) = 0.25 * ((1.0 + r) * (1.0 - s) * (1.0 - pow(t, 2.0)));
		m(0, 18) = 0.25 * ((1.0 + r) * (1.0 + s) * (1.0 - pow(t, 2.0)));
		m(0, 19) = 0.25 * ((1.0 - r) * (1.0 + s) * (1.0 - pow(t, 2.0)));
	}

	return N;
}

std::vector<double> Hexa20_weight()
{
	switch (Hexahedron20GPCount) {
	case 8: {
		return std::vector<double> (8, 1.0);
	}
	case 14: {
		std::vector<double> w(14, 0.0);
		double WF_scale_1 = 0.335180055401662;
		double WF_scale_2 = 0.886426592797784;
		for (int i = 0; i < 8; i++) {
			w[i] = (i < 8) ? WF_scale_1 : WF_scale_2;
		}
		return w;
	}
	default:
		std::cerr << "Unknown number of Tatrahedron10 GP count\n";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Hexahedron20::_dN = Hexa20_dN();
std::vector<DenseMatrix> Hexahedron20::_N = Hexa20_N();
std::vector<double> Hexahedron20::_weighFactor = Hexa20_weight();

bool Hexahedron20::match(eslocal *indices, eslocal n) {

#ifndef D3
	// Hexahedron20 is 3D element
	return false;
#endif

	if (n != 20) {
		return false;
	}

	for (eslocal i = 0; i < 19; i++) {
		for (eslocal j = i + 1; j < 20; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Hexahedron20::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	if (nodeIndex > 8) {
		result.resize(2);
	} else {
		result.resize(3);
	}

	switch (nodeIndex) {
	case 0: {
		result[0] = _indices[11];
		result[1] = _indices[8];
		result[2] = _indices[16];
		return result;
	}
	case 1: {
		result[0] = _indices[8];
		result[1] = _indices[9];
		result[2] = _indices[17];
		return result;
	}
	case 2: {
		result[0] = _indices[9];
		result[1] = _indices[10];
		result[2] = _indices[18];
		return result;
	}
	case 3: {
		result[0] = _indices[10];
		result[1] = _indices[11];
		result[2] = _indices[19];
		return result;
	}
	case 4: {
		result[0] = _indices[15];
		result[1] = _indices[12];
		result[2] = _indices[16];
		return result;
	}
	case 5: {
		result[0] = _indices[12];
		result[1] = _indices[13];
		result[2] = _indices[17];
		return result;
	}
	case 6: {
		result[0] = _indices[13];
		result[1] = _indices[14];
		result[2] = _indices[18];
		return result;
	}
	case 7: {
		result[0] = _indices[14];
		result[1] = _indices[15];
		result[2] = _indices[19];
		return result;
	}
	case 8: {
		result[0] = _indices[0];
		result[1] = _indices[1];
		return result;
	}
	case 9: {
		result[0] = _indices[1];
		result[1] = _indices[2];
		return result;
	}
	case 10: {
		result[0] = _indices[2];
		result[1] = _indices[3];
		return result;
	}
	case 11: {
		result[0] = _indices[3];
		result[1] = _indices[0];
		return result;
	}
	case 12: {
		result[0] = _indices[4];
		result[1] = _indices[5];
		return result;
	}
	case 13: {
		result[0] = _indices[5];
		result[1] = _indices[6];
		return result;
	}
	case 14: {
		result[0] = _indices[6];
		result[1] = _indices[7];
		return result;
	}
	case 15: {
		result[0] = _indices[7];
		result[1] = _indices[4];
		return result;
	}
	case 16: {
		result[0] = _indices[0];
		result[1] = _indices[4];
		return result;
	}
	case 17: {
		result[0] = _indices[1];
		result[1] = _indices[5];
		return result;
	}
	case 18: {
		result[0] = _indices[1];
		result[1] = _indices[6];
		return result;
	}
	case 19: {
		result[0] = _indices[3];
		result[1] = _indices[7];
		return result;
	}
	}
	return result;
}

std::vector<eslocal> Hexahedron20::getFace(size_t face) const
{
	std::vector<eslocal> result(4);

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

Hexahedron20::Hexahedron20(eslocal *indices)
{
	memcpy(_indices, indices, Hexahedron20NodesCount * sizeof(eslocal));
}


