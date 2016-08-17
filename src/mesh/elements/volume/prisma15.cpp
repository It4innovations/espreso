
#include "prisma15.h"
#include "../line/line3.h"
#include "../plane/triangle6.h"
#include "../plane/square8.h"

using namespace espreso;

size_t Prisma15::_counter = 0;

std::vector<Property> Prisma15::_DOFElement;
std::vector<Property> Prisma15::_DOFFace;
std::vector<Property> Prisma15::_DOFEdge;
std::vector<Property> Prisma15::_DOFPoint;
std::vector<Property> Prisma15::_DOFMidPoint;

static std::vector<std::vector< double> > Prisma15_rst()
{
	std::vector< std::vector<double> > rst(3, std::vector<double>(Prisma15GPCount));

	switch (Prisma15GPCount) {
	case 9: {
		double v1 = 1.0 / 6.0;
		double v2 = 4.0 / 6.0;
		double v3 = sqrt(3.0 / 5.0);
		double v4 = 0.0;
		rst[0] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
		rst[1] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
		rst[2] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };
		return rst;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Prisma15 GP count.";
		exit(EXIT_FAILURE);
	}
}


static std::vector<DenseMatrix> Prisma15_dN() {
	std::vector<DenseMatrix> dN(
		Prisma15GPCount,
		DenseMatrix(Point::size(), Prisma15NodesCount)
	);

	std::vector<std::vector< double> > rst = Prisma15_rst();

	for (unsigned int i = 0; i < Prisma15GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// dNr - derivation of basis function
		m(0, 0) = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
		m(0, 1) = ((t - 1.0) * (t - 2.0 * r + 2.0)) / 2.0 - r * (t - 1.0);
		m(0, 2) = 0.0;
		m(0, 3) = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
		m(0, 4) = r * (t + 1.0) + ((t + 1.0) * (2.0 * r + t - 2.0)) / 2.0;
		m(0, 5) = 0.0;
		m(0, 6) = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * r * (t - 1.0);
		m(0, 7) = (-2.0) * s * (t - 1.0);
		m(0, 8) = 2.0 * s * (t - 1.0);
		m(0, 9) = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * r * (t + 1.0);
		m(0, 10) = 2.0 * s * (t + 1.0);
		m(0, 11) = -2.0 * s * (t + 1.0);
		m(0, 12) = pow(t, 2.0) - 1.0;
		m(0, 13) = 1.0 - pow(t, 2.0);
		m(0, 14) = 0.0;

		// dNs - derivation of basis function
		m(1, 0) = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
		m(1, 1) = 0.0;
		m(1, 2) = ((t - 1.0) * (t - 2.0 * s + 2.0)) / 2.0 - s * (t - 1.0);
		m(1, 3) = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
		m(1, 4) = 0.0;
		m(1, 5) = s * (t + 1.0) + ((t + 1.0) * (2.0 * s + t - 2.0)) / 2.0;
		m(1, 6) = 2.0 * r * (t - 1.0);
		m(1, 7) = (-2.0) * r * (t - 1.0);
		m(1, 8) = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * s * (t - 1.0);
		m(1, 9) = (-2.0) * r * (t + 1.0);
		m(1, 10) = 2.0 * r * (t + 1.0);
		m(1, 11) = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * s * (t + 1.0);
		m(1, 12) = pow(t, 2.0) - 1.0;
		m(1, 13) = 0.0;
		m(1, 14) = 1.0 - pow(t, 2.0);

		// dNt - derivation of basis function
		m(2, 0) = -((r + s - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0 - ((t - 1.0) * (r + s - 1.0)) / 2.0;
		m(2, 1) = (r * (t - 2.0 * r + 2.0)) / 2.0 + (r * (t - 1.0)) / 2.0;
		m(2, 2) = (s * (t - 2.0 * s + 2.0)) / 2.0 + (s * (t - 1.0)) / 2.0;
		m(2, 3) = ((r + s - 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 - ((t + 1.0) * (r + s - 1.0)) / 2.0;
		m(2, 4) = (r * (2.0 * r + t - 2.0)) / 2.0 + (r * (t + 1.0)) / 2.0;
		m(2, 5) = (s * (2.0 * s + t - 2.0)) / 2.0 + (s * (t + 1.0)) / 2.0;
		m(2, 6) = 2.0 * r * (r + s - 1.0);
		m(2, 7) = (-2.0) * r * s;
		m(2, 8) = 2.0 * s * (r + s - 1.0);
		m(2, 9) = (-2.0) * r * (r + s - 1.0);
		m(2, 10) = 2.0 * r * s;
		m(2, 11) = (-2.0) * s * (r + s - 1.0);
		m(2, 12) = 2.0 * t * (r + s - 1.0);
		m(2, 13) = (-2.0) * r * t;
		m(2, 14) = (-2.0) * s * t;
	}

	return dN;
}

static std::vector<DenseMatrix> Prisma15_N() {
	std::vector<DenseMatrix> N(
		Prisma15GPCount,
		DenseMatrix(1, Prisma15NodesCount)
	);

	std::vector<std::vector< double> > rst = Prisma15_rst();

	for (unsigned int i = 0; i < Prisma15GPCount; i++) {
		DenseMatrix &m = N[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// basis function
		m(0, 0) = -(1.0 - r - s) * (1.0 - t) * (2.0 * r + 2.0 * s + t) / 2.0;
		m(0, 1) = r * (1.0 - t) * (2.0 * r - t - 2.0) / 2.0;
		m(0, 2) = s * (1.0 - t) * (2.0 * s - t - 2.0) / 2.0;
		m(0, 3) = -(1.0 - r - s) * (1.0 + t) * (2.0 * r + 2.0 * s - t) / 2.0;
		m(0, 4) = r * (t + 1.0) * (2.0 * r + t - 2.0) / 2.0;
		m(0, 5) = s * (t + 1.0) * (2.0 * s + t - 2.0) / 2.0;
		m(0, 6) = 2.0 * r * (1.0 - r - s) * (1.0 - t);
		m(0, 7) = 2.0 * r * s * (1.0 - t);
		m(0, 8) = 2.0 * s * (1.0 - r - s) * (1.0 - t);
		m(0, 9) = 2.0 * r * (1.0 - r - s) * (1.0 + t);
		m(0, 10) = 2.0 * r * s * (1.0 + t);
		m(0, 11) = 2.0 * s * (1.0 - r - s) * (1.0 + t);
		m(0, 12) = (1.0 - r - s) * (1.0 - pow(t, 2.0));
		m(0, 13) = r * (1.0 - pow(t, 2.0));
		m(0, 14) = s * (1.0 - pow(t, 2.0));
	}

	return N;
}

static std::vector<double> Prisma15_weight()
{
	std::vector<double> w;
	switch (Prisma15GPCount) {
	case 9: {
		double v1 = 5.0 / 54.0;
		double v2 = 8.0 / 54.0;
		w = { v1, v1, v1, v2, v2, v2, v1, v1, v1 };
		break;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Prisma15 GP count.";
		exit(EXIT_FAILURE);
	}
	return w;
}

std::vector<DenseMatrix> Prisma15::_dN = Prisma15_dN();
std::vector<DenseMatrix> Prisma15::_N = Prisma15_N();
std::vector<double> Prisma15::_weighFactor = Prisma15_weight();

bool Prisma15::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Prisma15 is 3D element
	return false;
#endif

	switch (n) {
	case Prisma15NodesCount:
		for (eslocal i = 0; i < Prisma15NodesCount - 1; i++) {
			for (eslocal j = i + 1; j < Prisma15NodesCount; j++) {
				if (Element::match(indices, i, j)) {
					return false;
				}
			}
		}
		return true;
	case 20:
		if (!Element::match(indices, 2, 3)) {
			return false;
		}
		if (!Element::match(indices, 3, 10)) {
			return false;
		}
		if (!Element::match(indices, 6, 7)) {
			return false;
		}
		if (!Element::match(indices, 7, 14)) {
			return false;
		}
		if (!Element::match(indices, 18, 19)) {
			return false;
		}

		eslocal various[Prisma15NodesCount] = { 0, 1, 2, 4, 5, 6, 8, 9, 11, 12, 13, 15, 16, 17, 18 };
		for (eslocal i = 0; i < Prisma15NodesCount - 1; i++) {
			for (eslocal j = i + 1; j < Prisma15NodesCount; j++) {
				if (Element::match(indices, various[i], various[j])) {
					return false;
				}
			}
		}
		return true;
	default:
		return false;
	}
}

std::vector<eslocal> Prisma15::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	if (nodeIndex > 5) {
		result.resize(2);
	} else {
		result.resize(3);
	}

	switch (nodeIndex) {
	case 0: {
		result[0] = _indices[6];
		result[1] = _indices[8];
		result[2] = _indices[12];
		return result;
	}
	case 1: {
		result[0] = _indices[6];
		result[1] = _indices[7];
		result[2] = _indices[13];
		return result;
	}
	case 2: {
		result[0] = _indices[7];
		result[1] = _indices[8];
		result[2] = _indices[14];
		return result;
	}
	case 3: {
		result[0] = _indices[9];
		result[1] = _indices[11];
		result[2] = _indices[12];
		return result;
	}
	case 4: {
		result[0] = _indices[9];
		result[1] = _indices[10];
		result[2] = _indices[13];
		return result;
	}
	case 5: {
		result[0] = _indices[10];
		result[1] = _indices[11];
		result[2] = _indices[14];
		return result;
	}
	case 6: {
		result[0] = _indices[0];
		result[1] = _indices[1];
		return result;
	}
	case 7: {
		result[0] = _indices[1];
		result[1] = _indices[2];
		return result;
	}
	case 8: {
		result[0] = _indices[0];
		result[1] = _indices[2];
		return result;
	}
	case 9: {
		result[0] = _indices[3];
		result[1] = _indices[4];
		return result;
	}
	case 10: {
		result[0] = _indices[4];
		result[1] = _indices[5];
		return result;
	}
	case 11: {
		result[0] = _indices[3];
		result[1] = _indices[5];
		return result;
	}
	case 12: {
		result[0] = _indices[0];
		result[1] = _indices[3];
		return result;
	}
	case 13: {
		result[0] = _indices[1];
		result[1] = _indices[4];
		return result;
	}
	case 14: {
		result[0] = _indices[2];
		result[1] = _indices[5];
		return result;
	}
	}
	return result;
}

void Prisma15::fillEdges()
{
	eslocal line[Line3NodesCount];
	_edges.reserve(Prisma15EdgeCount);

	for (size_t edge = 0; edge < 3; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 3];
		line[2] = _indices[ edge + 6     ];
		_edges.push_back(new Line3(line));
		_edges.back()->parentElements().push_back(this);

		line[0] = _indices[ edge          +  3];
		line[1] = _indices[(edge + 1) % 3 +  3];
		line[2] = _indices[ edge          + 9];
		_edges.push_back(new Line3(line));
		_edges.back()->parentElements().push_back(this);

		line[0] = _indices[edge     ];
		line[1] = _indices[edge +  3];
		line[2] = _indices[edge + 12];
		_edges.push_back(new Line3(line));
		_edges.back()->parentElements().push_back(this);
	}
}

void Prisma15::fillFaces()
{
	eslocal square[Square8NodesCount];
	eslocal triangle[Triangle6NodesCount];
	_faces.reserve(Prisma15FacesCount);

	for (size_t face = 0; face < 3; face++) {
		square[0] = _indices[ face              ];
		square[1] = _indices[(face + 1) % 3     ];
		square[2] = _indices[(face + 1) % 3 + 3 ];
		square[3] = _indices[ face          + 3 ];

		square[4] = _indices[ face          + 6 ];
		square[5] = _indices[(face + 1) % 3 + 12];
		square[6] = _indices[ face          + 9 ];
		square[7] = _indices[ face          + 12];
		_faces.push_back(new Square8(square));
		_faces.back()->parentElements().push_back(this);
	}

	triangle[0] = _indices[1];
	triangle[1] = _indices[0];
	triangle[2] = _indices[2];

	triangle[3] = _indices[6];
	triangle[4] = _indices[8];
	triangle[5] = _indices[7];
	_faces.push_back(new Triangle6(triangle));
	_faces.back()->parentElements().push_back(this);

	triangle[0] = _indices[3];
	triangle[1] = _indices[4];
	triangle[2] = _indices[5];

	triangle[3] = _indices[9];
	triangle[4] = _indices[10];
	triangle[5] = _indices[11];
	_faces.push_back(new Triangle6(triangle));
	_faces.back()->parentElements().push_back(this);
}

Prisma15::Prisma15(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 20:
		_indices[0] = indices[0];
		_indices[1] = indices[1];
		_indices[2] = indices[2];
		_indices[3] = indices[4];
		_indices[4] = indices[5];
		_indices[5] = indices[6];
		_indices[6] = indices[8];
		_indices[7] = indices[9];
		_indices[8] = indices[11];
		_indices[9] = indices[12];
		_indices[10] = indices[13];
		_indices[11] = indices[15];
		_indices[12] = indices[16];
		_indices[13] = indices[17];
		_indices[14] = indices[19];
		break;
	case 15:
		memcpy(_indices, indices, 15 * sizeof(eslocal));
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Prisma15 from " << n << " elements.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Prisma15::Prisma15(std::ifstream &is)
{
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}





