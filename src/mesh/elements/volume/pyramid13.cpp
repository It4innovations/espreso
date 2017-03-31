
#include <cstring>
#include <fstream>

#include "pyramid13.h"
#include "../line/line3.h"
#include "../plane/triangle6.h"
#include "../plane/square8.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

std::vector<Property> Pyramid13::_DOFElement;
std::vector<Property> Pyramid13::_DOFFace;
std::vector<Property> Pyramid13::_DOFEdge;
std::vector<Property> Pyramid13::_DOFPoint;
std::vector<Property> Pyramid13::_DOFMidPoint;

std::vector<std::vector<eslocal> > Pyramid13::_facesNodes = {
	{ 3, 2, 1, 0,  7,  6, 5, 8 },
	{ 0, 1, 4, 5, 10,  9 },
	{ 1, 2, 4, 6, 11, 10 },
	{ 2, 3, 4, 7, 12, 11 },
	{ 3, 0, 4, 8,  9, 12 }
};
std::vector<std::vector<eslocal> > Pyramid13::_edgesNodes = {
	{ 0, 1,  5 },
	{ 1, 2,  6 },
	{ 2, 3,  7 },
	{ 3, 0,  8 },
	{ 0, 4,  9 },
	{ 1, 4, 10 },
	{ 2, 4, 11 },
	{ 3, 4, 12 }
};

static std::vector< std::vector< double> > Pyramid13_rst()
{
	std::vector< std::vector<double> > rst(3, std::vector<double>(Pyramid13GPCount));

	switch (Pyramid13GPCount) {
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
		ESINFO(ERROR) << "Unknown number of Pyramid13 GP count.";
		exit(EXIT_FAILURE);
	}
}

static std::vector<DenseMatrix> Pyramid13_dN() {
	std::vector<DenseMatrix> dN(
		Pyramid13GPCount,
		DenseMatrix(3, Pyramid13NodesCount)
	);

	std::vector< std::vector< double> > _pyramid13_rst = Pyramid13_rst();

	for (unsigned int i = 0; i < Pyramid13GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		double r = _pyramid13_rst[0][i];
		double s = _pyramid13_rst[1][i];
		double t = _pyramid13_rst[2][i];

		// dNr - derivation of basis function
		m(0, 0)  = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
		m(0, 1)  = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
		m(0, 2)  =  (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0) + (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
		m(0, 3)  =  (t / 2.0 - 0.5) * (s + 1.0) * (r - 1.0) * (t / 8.0 - 1.0 / 8.0) - (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0);
		m(0, 4)  =  0.0;
		m(0, 5)  =  r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s - 1.0);
		m(0, 6)  = -((pow(s, 2.0) - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
		m(0, 7)  = -r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s + 1.0);
		m(0, 8)  =  ((pow(s, 2.0) - 1) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
		m(0, 9)  = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
		m(0, 10) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
		m(0, 11) = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);
		m(0, 12) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);
		//  dNs - derivation of basis function
		m(1, 0)  = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) - 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
		m(1, 1)  =  (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0) - (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) - s * (t / 2.0 - 1.0 / 2.0) + 1.0);
		m(1, 2)  =  (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) + 1.0) + (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
		m(1, 3)  = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s * (t / 2.0 - 1.0 / 2.0) - r * (t / 2.0 - 1.0 / 2.0) + 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0);
		m(1, 4)  =  0.0;
		m(1, 5)  =  ((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
		m(1, 6)  = -s * (t / 2.0 - 1.0 / 2.0) * (t / 2.0 - 0.5) * (r + 1.0);
		m(1, 7)  = -((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
		m(1, 8)  =  s * (t / 2.0 - 0.5) * (t / 2.0 - 0.5) * (r - 1.0);
		m(1, 9)  = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);
		m(1, 10) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
		m(1, 11) = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
		m(1, 12) =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);
		//  dNt - derivation of basis function
		m(2, 0)  = -((r - 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0)) / 8.0 - (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
		m(2, 1)  = -((r + 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0)) / 8.0 - (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
		m(2, 2)  =  ((r + 1.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0)) / 8.0 + (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
		m(2, 3)  =  (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0) - ((r - 1.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0)) / 8.0;
		m(2, 4)  =  t + 0.5;
		m(2, 5)  =  ((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (s - 1.0)) / 2.0;
		m(2, 6)  = -((pow(s, 2.0) - 1.0) * (t / 2.0 - 0.5) * (r + 1.0)) / 2.0;
		m(2, 7)  = -((pow(r, 2.0) - 1.0) * (t / 2.0 - 0.5) * (s + 1.0)) / 2.0;
		m(2, 8)  =  ((pow(s, 2.0) - 1.0) * (t / 2.0 - 0.5) * (r - 1.0)) / 2.0;
		m(2, 9)  =  ((t / 2.0 - 0.5) * (r + s - r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r + s - r * s - 1.0)) / 2.0;
		m(2, 10) = -((t / 2.0 - 0.5) * (r - s - r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r - s - r * s + 1.0)) / 2.0;
		m(2, 11) = -((t / 2.0 - 0.5) * (r + s + r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r + s + r * s + 1.0)) / 2.0;
		m(2, 12) =  ((t / 2.0 - 0.5) * (r - s + r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r - s + r * s - 1.0)) / 2.0;
	}

	return dN;
}

static std::vector<DenseMatrix> Pyramid13_N() {
	std::vector<DenseMatrix> N(
		Pyramid13GPCount,
		DenseMatrix(1, Pyramid13NodesCount)
	);

	std::vector< std::vector< double> > _pyramid13_rst = Pyramid13_rst();

	for (unsigned int i = 0; i < Pyramid13GPCount; i++) {
		DenseMatrix &m = N[i];

		double r = _pyramid13_rst[0][i];
		double s = _pyramid13_rst[1][i];
		double t = _pyramid13_rst[2][i];

		// basis function
		m(0, 0)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 - s) * (-1.0 - (0.5 * (1.0 - t)) * r - (0.5 * (1.0 - t)) * s));
		m(0, 1)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 - s) * (-1.0 + (0.5 * (1 - t))   * r - (0.5 * (1.0 - t)) * s));
		m(0, 2)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 + s) * (-1.0 + (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
		m(0, 3)  = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 + s) * (-1.0 - (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
		m(0, 4)  = (1.0 - (0.5 * (1.0 - t))) * (1.0 - 2.0 * (0.5 * (1.0 - t)));
		m(0, 5)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - s) * (1.0 - pow(r, 2.0));
		m(0, 6)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + r) * (1.0 - pow(s, 2.0));
		m(0, 7)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + s) * (1.0 - pow(r, 2.0));
		m(0, 8)  = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - r) * (1.0 - pow(s, 2.0));
		m(0, 9)  = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r - s + r * s);
		m(0, 10) = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r - s - r * s);
		m(0, 11) = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r + s + r * s);
		m(0, 12) = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r + s - r * s);
	}
	return N;
}

static std::vector<double> Pyramid13_weight()
{
	switch (Pyramid13GPCount) {
	case 8: {
		return std::vector<double> (8, 1.0);
	}
	case 14: {
		std::vector<double> w(8, 0.335180055401662);
		w.resize(14, 0.886426592797784);
		return w;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Pyramid13 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Pyramid13::_dN = Pyramid13_dN();
std::vector<DenseMatrix> Pyramid13::_N = Pyramid13_N();
std::vector<double> Pyramid13::_weighFactor = Pyramid13_weight();

bool Pyramid13::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Hexahedron8 is 3D element
	return false;
#endif

	switch (n) {
	case 13:
		for (eslocal i = 0; i < 12; i++) {
			for (eslocal j = i + 1; j < 13; j++) {
				if (Element::match(indices, i, j)) {
					return false;
				}
			}
		}
		return true;
	case 20: {
		if (!Element::match(indices, 4, 5)) {
			return false;
		}
		if (!Element::match(indices, 5, 6)) {
			return false;
		}
		if (!Element::match(indices, 6, 7)) {
			return false;
		}

		eslocal various[5] = { 0, 1, 2, 3, 4};
		for (eslocal i = 0; i < 4; i++) {
			for (eslocal j = i + 1; j < 5; j++) {
				if (Element::match(indices, various[i], various[j])) {
					return false;
				}
			}
		}
		return true;
	}
	default:
		return false;
	}

}

std::vector<eslocal> Pyramid13::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	if (nodeIndex < 4) {
		result.resize(3);
	} else if (nodeIndex>4) {
		result.resize(2);
	}else {
    result.resize(4);
  }


	switch (nodeIndex) {
	case 0: {
		result[0] = _indices[5];
		result[1] = _indices[8];
		result[2] = _indices[9];
		return result;
	}
	case 1: {
		result[0] = _indices[5];
		result[1] = _indices[6];
		result[2] = _indices[10];
		return result;
	}
	case 2: {
		result[0] = _indices[6];
		result[1] = _indices[7];
		result[2] = _indices[11];
		return result;
	}
	case 3: {
		result[0] = _indices[7];
		result[1] = _indices[8];
		result[2] = _indices[12];
		return result;
	}
	case 4: {
		result[0] = _indices[9];
		result[1] = _indices[10];
		result[2] = _indices[11];
		result[3] = _indices[12];
		return result;
	}
	case 5: {
		result[0] = _indices[0];
		result[1] = _indices[1];
		return result;
	}
	case 6: {
		result[0] = _indices[1];
		result[1] = _indices[2];
		return result;
	}
	case 7: {
		result[0] = _indices[2];
		result[1] = _indices[3];
		return result;
	}
	case 8: {
		result[0] = _indices[0];
		result[1] = _indices[3];
		return result;
	}
	case 9: {
		result[0] = _indices[0];
		result[1] = _indices[4];
		return result;
	}
	case 10: {
		result[0] = _indices[1];
		result[1] = _indices[4];
		return result;
	}
	case 11: {
		result[0] = _indices[2];
		result[1] = _indices[4];
		return result;
	}
	case 12: {
		result[0] = _indices[3];
		result[1] = _indices[4];
		return result;
	}
	}
	return result;
}

size_t Pyramid13::fillEdges()
{
	eslocal line[Line3NodesCount];

	if (_edges.size() == Pyramid13EdgeCount) {
		return Pyramid13EdgeCount;
	}
	_edges.reserve(Pyramid13EdgeCount);

	size_t filled = _edges.size();

	for (size_t e = 0 ; e < Pyramid13EdgeCount; e++) {
		for (size_t n = 0; n < Line3NodesCount; n++) {
			line[n] = _indices[_edgesNodes[e][n]];
		}
		addUniqueEdge<Line3>(line, filled, Line2NodesCount);
	}

	return filled;
}

size_t Pyramid13::fillFaces()
{
	eslocal square[Square8NodesCount];
	eslocal triangle[Triangle6NodesCount];

	if (_faces.size() == Pyramid13FacesCount) {
		return Pyramid13FacesCount;
	}
	_faces.reserve(Pyramid13FacesCount);

	size_t filled = _faces.size();

	for (size_t n = 0; n < Square8NodesCount; n++) {
		square[n] = _indices[_facesNodes[0][n]];
	}
	addUniqueFace<Square8>(square, filled, Square4NodesCount);

	for (size_t f = 1 ; f < Pyramid13FacesCount; f++) {
		for (size_t n = 0; n < Triangle6NodesCount; n++) {
			triangle[n] = _indices[_facesNodes[f][n]];
		}
		addUniqueFace<Triangle6>(triangle, filled, Triangle3NodesCount);
	}

	return filled;
}

Element* Pyramid13::addFace(const std::vector<eslocal> &nodes)
{
	for (size_t f = 0; f < faces(); f++) {
		size_t found;
		for (found = 0; found < _facesNodes[f].size(); found++) {
			if (!std::binary_search(nodes.begin(), nodes.end(), _indices[_facesNodes[f][found]])) {
				break;
			}
		}
		if (found == _facesNodes[f].size()) {
			if (f == 0) {
				eslocal square[Square8NodesCount];
				for (size_t n = 0; n < Square8NodesCount; n++) {
					square[n] = _indices[_facesNodes[f][n]];
				}
				return addUniqueFace<Square8>(square, _faces.size(), Square4NodesCount);
			} else {
				eslocal triangle[Triangle6NodesCount];
				for (size_t n = 0; n < Triangle6NodesCount; n++) {
					triangle[n] = _indices[_facesNodes[f][n]];
				}
				return addUniqueFace<Triangle6>(triangle, _faces.size(), Triangle3NodesCount);
			}
		}
	}
	return NULL;
}

Pyramid13::Pyramid13(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 13:
		memcpy(_indices, indices, 13 * sizeof(eslocal));
		break;
	case 20:
		_indices[ 0] = indices[ 0];
		_indices[ 1] = indices[ 1];
		_indices[ 2] = indices[ 2];
		_indices[ 3] = indices[ 3];
		_indices[ 4] = indices[ 4];
		_indices[ 5] = indices[ 8];
		_indices[ 6] = indices[ 9];
		_indices[ 7] = indices[10];
		_indices[ 8] = indices[11];
		_indices[ 9] = indices[16];
		_indices[10] = indices[17];
		_indices[11] = indices[18];
		_indices[12] = indices[19];
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Pyramid13 from " << n << " indices.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Pyramid13::Pyramid13(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}






