
#include <cstring>
#include <fstream>

#include "hexahedron20.h"
#include "../line/line3.h"
#include "../plane/square8.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

std::vector<Property> Hexahedron20::_DOFElement;
std::vector<Property> Hexahedron20::_DOFFace;
std::vector<Property> Hexahedron20::_DOFEdge;
std::vector<Property> Hexahedron20::_DOFPoint;
std::vector<Property> Hexahedron20::_DOFMidPoint;

std::vector<std::vector<eslocal> > Hexahedron20::_facesNodes = {
	{ 0, 1, 5, 4,  8, 17, 12, 16 },
	{ 3, 2, 1, 0, 10,  9,  8, 11 },
	{ 4, 5, 6, 7, 12, 13, 14, 15 },
	{ 7, 6, 2, 3, 14, 18, 10, 19 },
	{ 1, 2, 6, 5,  9, 18, 13, 17 },
	{ 3, 0, 4, 7, 11, 16, 15, 19 }
};

std::vector<std::vector<eslocal> > Hexahedron20::_edgesNodes = {
	{ 0, 1,  8 },
	{ 1, 2,  9 },
	{ 2, 3, 10 },
	{ 3, 0, 11 },
	{ 4, 5, 12 },
	{ 5, 6, 13 },
	{ 6, 7, 14 },
	{ 7, 4, 15 },
	{ 0, 4, 16 },
	{ 1, 5, 17 },
	{ 2, 6, 18 },
	{ 3, 7, 19 }
};

static std::vector<std::vector< double> > Hexa20_rst()
{
	std::vector< std::vector<double> > rst(3, std::vector<double>(Hexahedron20GPCount));

	switch (Hexahedron20GPCount) {
	case 8: {
		double v = 0.577350269189625953;
		rst[0] = {  v,  v,  v,  v, -v, -v, -v, -v };
		rst[1] = { -v, -v,  v,  v, -v, -v,  v,  v };
		rst[2] = { -v,  v, -v,  v, -v,  v, -v,  v };
		return rst;
	}
	case 14: {
		double v1 = 0.758786910639329015;
		double v2 = 0.795822425754222018;
		double v3 = 0;
		rst[0] = { -v1,  v1,  v1, -v1, -v1,  v1,  v1, -v1,  v3,  v3,  v2, v3, -v2, v3 };
		rst[1] = { -v1, -v1,  v1,  v1, -v1, -v1,  v1,  v1,  v3, -v2,  v3, v2,  v3, v3 };
		rst[2] = { -v1, -v1, -v1, -v1,  v1,  v1,  v1,  v1, -v2,  v3,  v3, v3,  v3, v2 };
		return rst;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Hexahedron20 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Hexa20_dN()
{
	std::vector<DenseMatrix> dN(
		Hexahedron20GPCount,
		DenseMatrix(3, Hexahedron20NodesCount)
	);

	std::vector<std::vector< double> > rst = Hexa20_rst();

	for (unsigned int i = 0; i < Hexahedron20GPCount; i++) {
		DenseMatrix &m = dN[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];


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

static std::vector<DenseMatrix> Hexa20_N() {
	std::vector<DenseMatrix> N(
		Hexahedron20GPCount,
		DenseMatrix(1, Hexahedron20NodesCount)
	);

	std::vector< std::vector< double> > rst = Hexa20_rst();

	for (unsigned int i = 0; i < Hexahedron20GPCount; i++) {
		DenseMatrix &m = N[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

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

static std::vector<double> Hexa20_weight()
{
	switch (Hexahedron20GPCount) {
	case 8: {
		return std::vector<double> (8, 1.0);
	}
	case 14: {
		std::vector<double> w(8, 0.335180055401662);
		w.resize(14, 0.886426592797784);
		return w;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Hexahedron20 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Hexahedron20::_dN = Hexa20_dN();
std::vector<DenseMatrix> Hexahedron20::_N = Hexa20_N();
std::vector<double> Hexahedron20::_weighFactor = Hexa20_weight();

bool Hexahedron20::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Hexahedron20 is 3D element
	return false;
#endif

	switch (n) {
	case Hexahedron20NodesCount:
		for (eslocal i = 0; i < Hexahedron20NodesCount - 1; i++) {
			for (eslocal j = i + 1; j < Hexahedron20NodesCount; j++) {
				if (Element::match(indices, i, j)) {
					return false;
				}
			}
		}
		return true;
	default:
		return false;
	}
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

size_t Hexahedron20::fillEdges()
{
	eslocal line[Line3NodesCount];
	if (_edges.size() == Hexahedron20EdgeCount) {
		return Hexahedron20EdgeCount;
	}

	_edges.reserve(Hexahedron20EdgeCount);

	size_t filled = _edges.size();

	for (size_t e = 0 ; e < Hexahedron20EdgeCount; e++) {
		for (size_t n = 0; n < Line3NodesCount; n++) {
			line[n] = _indices[_edgesNodes[e][n]];
		}
		addUniqueEdge<Line3>(line, filled, Line2NodesCount);
	}

	return filled;
}

size_t Hexahedron20::fillFaces()
{
	eslocal square[Square8NodesCount];
	if (_faces.size() == Hexahedron20FacesCount) {
		return Hexahedron20FacesCount;
	}
	_faces.reserve(Hexahedron20FacesCount);

	size_t filled = _faces.size();

	for (size_t f = 0 ; f < Hexahedron20FacesCount; f++) {
		for (size_t n = 0; n < Square8NodesCount; n++) {
			square[n] = _indices[_facesNodes[f][n]];
		}
		addUniqueFace<Square8>(square, filled, Square4NodesCount);
	}
	return filled;
}

Element* Hexahedron20::addFace(const std::vector<eslocal> &nodes)
{
	for (size_t f = 0; f < faces(); f++) {
		size_t found;
		for (found = 0; found < _facesNodes[f].size(); found++) {
			if (!std::binary_search(nodes.begin(), nodes.end(), _indices[_facesNodes[f][found]])) {
				break;
			}
		}
		if (found == _facesNodes[f].size()) {
			eslocal square[Square8NodesCount];
			for (size_t n = 0; n < Square8NodesCount; n++) {
				square[n] = _indices[_facesNodes[f][n]];
			}
			return addUniqueFace<Square8>(square, _faces.size(), Square4NodesCount);
		}
	}
	return NULL;
}

Hexahedron20::Hexahedron20(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 20:
		memcpy(_indices, indices, Hexahedron20NodesCount * sizeof(eslocal));
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Hexahedron20 from " << n << " indices.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Hexahedron20::Hexahedron20(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}


