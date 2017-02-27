
#include <cstring>
#include <fstream>

#include "hexahedron8.h"
#include "../line/line2.h"
#include "../plane/square4.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

std::vector<Property> Hexahedron8::_DOFElement;
std::vector<Property> Hexahedron8::_DOFFace;
std::vector<Property> Hexahedron8::_DOFEdge;
std::vector<Property> Hexahedron8::_DOFPoint;
std::vector<Property> Hexahedron8::_DOFMidPoint;

std::vector<std::vector<eslocal> > Hexahedron8::_facesNodes = {
	{ 0, 1, 5, 4 },
	{ 3, 2, 1, 0 },
	{ 4, 5, 6, 7 },
	{ 7, 6, 2, 3 },
	{ 1, 2, 6, 5 },
	{ 3, 0, 4, 7 }
};

std::vector<std::vector<eslocal> > Hexahedron8::_edgesNodes = {
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 3 },
	{ 3, 0 },
	{ 4, 5 },
	{ 5, 6 },
	{ 6, 7 },
	{ 7, 4 },
	{ 0, 4 },
	{ 1, 5 },
	{ 2, 6 },
	{ 3, 7 }
};

static std::vector<DenseMatrix> Hexa_dN() {
	std::vector<DenseMatrix> dN(
		Hexahedron8GPCount,
		DenseMatrix(3, Hexahedron8NodesCount)
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

static std::vector<DenseMatrix> Hexa_N() {
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
std::vector<double> Hexahedron8::_weighFactor(Hexahedron8GPCount, 1);

bool Hexahedron8::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Hexahedron8 is 3D element
	return false;
#endif

	switch (n) {
	case Hexahedron8NodesCount:
		for (eslocal i = 0; i < Hexahedron8NodesCount - 1; i++) {
			for (eslocal j = i + 1; j < Hexahedron8NodesCount; j++) {
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

std::vector<eslocal> Hexahedron8::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(3);

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

size_t Hexahedron8::fillEdges()
{
	eslocal line[Line2NodesCount];

	if (_edges.size() == Hexahedron8EdgeCount) {
		return Hexahedron8EdgeCount;
	}
	_edges.reserve(Hexahedron8EdgeCount);

	size_t filled = _edges.size();

	for (size_t e = 0 ; e < Hexahedron8EdgeCount; e++) {
		for (size_t n = 0; n < Line2NodesCount; n++) {
			line[n] = _indices[_edgesNodes[e][n]];
		}
		addUniqueEdge<Line2>(line, filled, Line2NodesCount);
	}

	return filled;
}

size_t Hexahedron8::fillFaces()
{
	eslocal square[Square4NodesCount];
	if (_faces.size() == Hexahedron8FacesCount) {
		return Hexahedron8FacesCount;
	}
	_faces.reserve(Hexahedron8FacesCount);

	size_t filled = _faces.size();

	for (size_t f = 0 ; f < Hexahedron8FacesCount; f++) {
		for (size_t n = 0; n < Square4NodesCount; n++) {
			square[n] = _indices[_facesNodes[f][n]];
		}
		addUniqueFace<Square4>(square, filled, Square4NodesCount);
	}

	return filled;
}

Hexahedron8::Hexahedron8(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 8:
		memcpy(_indices, indices, Hexahedron8NodesCount * sizeof(eslocal));
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Hexahedron8 from " << n << " indices.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Hexahedron8::Hexahedron8(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}


