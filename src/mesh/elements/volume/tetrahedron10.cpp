
#include <cstring>
#include <fstream>

#include "tetrahedron10.h"
#include "../line/line3.h"
#include "../plane/triangle6.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

std::vector<Property> Tetrahedron10::_DOFElement;
std::vector<Property> Tetrahedron10::_DOFFace;
std::vector<Property> Tetrahedron10::_DOFEdge;
std::vector<Property> Tetrahedron10::_DOFPoint;
std::vector<Property> Tetrahedron10::_DOFMidPoint;

std::vector<std::vector<eslocal> > Tetrahedron10::_facesNodes = {
	{ 0, 1, 3, 4, 8, 7 },
	{ 1, 2, 3, 5, 9, 8 },
	{ 2, 0, 3, 6, 7, 9 },
	{ 2, 1, 0, 5, 4, 6 }
};

std::vector<std::vector<eslocal> > Tetrahedron10::_edgesNodes = {
	{ 0, 1, 4 },
	{ 1, 2, 5 },
	{ 2, 0, 6 },
	{ 0, 3, 7 },
	{ 1, 3, 8 },
	{ 2, 3, 9 }
};

static std::vector<DenseMatrix> Tetra10_dN()
{
	std::vector<DenseMatrix> dN(
		Tetrahedron10GPCount,
		DenseMatrix(3, Tetrahedron10NodesCount)
	);

	std::vector<double> rv;
	std::vector<double> sv;
	std::vector<double> tv;

	// WARNINKG: use only 15 gausse points
	switch (Tetrahedron10GPCount) {
	case 4:
		rv = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		sv = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		tv = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
		break;
	case 5:
		rv = {0.2500000000000000, 0.5000000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
		sv = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000};
		tv = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000, 0.1666666666666667};
		break;
	case 11:
		rv = { 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0, 0.25 };
		sv = { 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.25 };
		tv = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.25 };
		break;
	case 15:
		rv = {
				0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
				0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
				0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
				0.0665501535736643, 0.4334498464263357, 0.4334498464263357};
		sv = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
				0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
				0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
				0.4334498464263357, 0.0665501535736643, 0.4334498464263357};
		tv = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
				0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
				0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
				0.4334498464263357, 0.4334498464263357, 0.0665501535736643};
		break;
	default:
		ESINFO(ERROR) << "Unknown number of Tatrahedron10 GP count.";
	}

	for (unsigned int i = 0; i < Tetrahedron10GPCount; i++) {
		double r = rv[i];
		double s = sv[i];
		double t = tv[i];

		DenseMatrix &m = dN[i];

		m(0, 0) = 4.0 * r - 1.0;
		m(0, 1) = 0;
		m(0, 2) = 0;
		m(0, 3) = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
		m(0, 4) = 4.0 * s;
		m(0, 5) = 0;
		m(0, 6) = 4.0 * t;
		m(0, 7) = -8.0 * r - 4.0 * s - 4.0 * t + 4.0;
		m(0, 8) = -4.0 * s;
		m(0, 9) = -4.0 * t;

		m(2, 0) = 0;
		m(2, 1) = 4.0 * s - 1.0;
		m(2, 2) = 0 ;
		m(2, 3) = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
		m(2, 4) = 4.0 * r;
		m(2, 5) = 4.0 * t;
		m(2, 6) = 0;
		m(2, 7) = -4.0 * r;
		m(2, 8) = -4.0 * r - 8.0 * s - 4.0 * t + 4.0;
		m(2, 9) = -4.0 * t;

		m(1, 0) = 0;
		m(1, 1) = 0;
		m(1, 2) = 4.0 * t - 1.0;
		m(1, 3) = 4.0 * r + 4.0 * s + 4.0* t  - 3.0;
		m(1, 4) = 0;
		m(1, 5) = 4.0 * s;
		m(1, 6) = 4.0 * r;
		m(1, 7) = -4.0 * r;
		m(1, 8) = -4.0 * s;
		m(1, 9) = -4.0 * r - 4.0 * s - 8.0 * t + 4.0;
	}

	return dN;
}

static std::vector<DenseMatrix> Tetra10_N() {
	std::vector<DenseMatrix> N(
		Tetrahedron10GPCount,
		DenseMatrix(1, Tetrahedron10NodesCount)
	);

	std::vector<double> rv;
	std::vector<double> sv;
	std::vector<double> tv;


	switch (Tetrahedron10GPCount) {
	case 4:
		rv = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		sv = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		tv = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
		break;
	case 5:
		rv = {0.2500000000000000, 0.5000000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
		sv = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000};
		tv = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000, 0.1666666666666667};
		break;
	case 11:
		rv = {
				0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714,
				0.0714285714285714, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992,
				0.3994035761667992, 0.1005964238332008, 0.1005964238332008};
		sv = {
				0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
				0.7857142857142857, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992,
				0.1005964238332008, 0.3994035761667992, 0.1005964238332008};
		tv = {
				0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857,
				0.0714285714285714, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008,
				0.1005964238332008, 0.1005964238332008, 0.3994035761667992};
		break;
	case 15:
		rv = {
				0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
				0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
				0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
				0.0665501535736643, 0.4334498464263357, 0.4334498464263357};
		sv = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
				0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
				0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
				0.4334498464263357, 0.0665501535736643, 0.4334498464263357};
		tv = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
				0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
				0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
				0.4334498464263357, 0.4334498464263357, 0.0665501535736643};
		break;
	default:
		ESINFO(ERROR) << "Unknown number of Tatrahedron10 GP count.";
	}

	for (unsigned int i = 0; i < Tetrahedron10GPCount; i++) {
		double r = rv[i];
		double s = sv[i];
		double t = tv[i];

		N[i](0, 0) = r * (2.0 * r - 1.0);
		N[i](0, 1) = s * (2.0 * s - 1.0);
		N[i](0, 2) = t * (2.0 * t - 1.0);
		N[i](0, 3) = 2.0 * r * r + 4.0 * r * s + 4.0 * r * t - 3.0 * r + 2.0* s * s + 4.0 * s * t - 3.0 * s + 2.0 * t * t - 3.0 * t + 1.0;
		N[i](0, 4) = 4.0 * r * s;
		N[i](0, 5) = 4.0 * s * t;
		N[i](0, 6) = 4.0 * r * t;
		N[i](0, 7) = r * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
		N[i](0, 8) = s * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
		N[i](0, 9) = t * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
	}

	return N;
}


static std::vector<double> Tetra10_Weight()
{
	switch (Tetrahedron10GPCount) {
	case 4: {
		return std::vector<double> (4, 1 / 24.0);
	}
	case 5: {
		std::vector<double> w(5, 3 / 40.0);
		w[0] = - 2 / 15.0;
		return w;
	}
	case 11: {
		std::vector<double> w(11);
		w[0] = w[1] = w[2] = w[3] = 0.002777778;
		w[5] = w[6] = w[7] = w[8] = w[9] = 0.011111111;
		w[10] = 0.0888888888889;
		return w;
	}
	case 15: {
		std::vector<double> w(15);
		w[0] = 0.030283678097089;
		w[1] = w[2] = w[3] = w[4] = 0.006026785714286;
		w[5] = w[6] = w[7] = w[8] = 0.011645249086029;
		w[9] = w[10] = w[11] = w[12] = w[13] = w[14] = 0.010949141561386;
		return w;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Tatrahedron10 GP count.";
		exit(EXIT_FAILURE);
	}
}



std::vector<DenseMatrix> Tetrahedron10::_dN = Tetra10_dN();
std::vector<DenseMatrix> Tetrahedron10::_N = Tetra10_N();
std::vector<double> Tetrahedron10::_weighFactor = Tetra10_Weight();

bool Tetrahedron10::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Tetrahedron10 is 3D element
	return false;
#endif

	switch (n) {
	case 20: {
		if (!Element::match(indices, 2, 3)) {
			return false;
		}
		if (!Element::match(indices, 2, 10)) {
			return false;
		}
		if (!Element::match(indices, 18, 19)) {
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

		eslocal various[10] = { 0, 1, 2, 4, 8, 9 ,11, 16, 17, 18 };
		for (eslocal i = 0; i < 9; i++) {
			for (eslocal j = i + 1; j < 10; j++) {
				if (Element::match(indices, various[i], various[j])) {
					return false;
				}
			}
		}
		return true;
	}
	case 10:
		for (eslocal i = 0; i < 9; i++) {
			for (eslocal j = i + 1; j < 10; j++) {
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

std::vector<eslocal> Tetrahedron10::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
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

size_t Tetrahedron10::fillEdges()
{
	eslocal line[Line3NodesCount];

	if (_edges.size() == Tetrahedron10EdgeCount) {
		return Tetrahedron10EdgeCount;
	}
	_edges.reserve(Tetrahedron10EdgeCount);

	size_t filled = _edges.size();

	for (size_t e = 0; e < Tetrahedron10EdgeCount; e++) {
		for (size_t n = 0; n < Line3NodesCount; n++) {
			line[n] = _indices[_edgesNodes[e][n]];
		}
		addUniqueEdge<Line3>(line, filled, Line2NodesCount);
	}

	return filled;
}

size_t Tetrahedron10::fillFaces()
{
	eslocal triangle[Triangle6NodesCount];

	if (_faces.size() == Tetrahedron10FacesCount) {
		return Tetrahedron10FacesCount;
	}
	_faces.reserve(Tetrahedron10FacesCount);

	size_t filled = _faces.size();

	for (size_t f = 0; f < Tetrahedron10FacesCount; f++) {
		for (size_t n = 0; n < Triangle6NodesCount; n++) {
			triangle[n] = _indices[_facesNodes[f][n]];
		}
		addUniqueFace<Triangle6>(triangle, filled, Triangle3NodesCount);
	}

	return filled;
}

Element* Tetrahedron10::addFace(const std::vector<eslocal> &nodes)
{
	for (size_t f = 0; f < faces(); f++) {
		size_t found;
		for (found = 0; found < _facesNodes[f].size(); found++) {
			if (!std::binary_search(nodes.begin(), nodes.end(), _indices[_facesNodes[f][found]])) {
				break;
			}
		}
		if (found == _facesNodes[f].size()) {
			eslocal triangle[Triangle6NodesCount];
			for (size_t n = 0; n < Triangle6NodesCount; n++) {
				triangle[n] = _indices[_facesNodes[f][n]];
			}
			return addUniqueFace<Triangle6>(triangle, _faces.size(), Triangle3NodesCount);
		}
	}
	return NULL;
}

Tetrahedron10::Tetrahedron10(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 10:
		memcpy(_indices, indices, 10 * sizeof(eslocal));
		break;
	case 20:
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
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Tetrahedron10 from " << n << " indices.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Tetrahedron10::Tetrahedron10(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}



