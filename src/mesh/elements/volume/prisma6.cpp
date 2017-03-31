
#include <cstring>
#include <fstream>

#include "prisma6.h"
#include "../line/line2.h"
#include "../plane/triangle3.h"
#include "../plane/square4.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

std::vector<Property> Prisma6::_DOFElement;
std::vector<Property> Prisma6::_DOFFace;
std::vector<Property> Prisma6::_DOFEdge;
std::vector<Property> Prisma6::_DOFPoint;
std::vector<Property> Prisma6::_DOFMidPoint;

std::vector<std::vector<eslocal> > Prisma6::_facesNodes = {
	{ 0, 1, 4, 3 },
	{ 1, 2, 5, 4 },
	{ 2, 0, 3, 5 },
	{ 2, 1, 0 },
	{ 3, 4, 5 }
};

std::vector<std::vector<eslocal> > Prisma6::_edgesNodes = {
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 0 },
	{ 3, 4 },
	{ 4, 5 },
	{ 5, 3 },
	{ 0, 3 },
	{ 1, 4 },
	{ 2, 5 }
};

static std::vector< std::vector< double> > Prisma6_rst()
{
	std::vector< std::vector<double> > rst(3, std::vector<double>(Prisma6GPCount));

	switch (Prisma6GPCount) {
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
		ESINFO(ERROR) << "Unknown number of Prisma6 GP count.";
		exit(EXIT_FAILURE);
	}
}

static std::vector<DenseMatrix> Prisma6_dN() {
	std::vector<DenseMatrix> dN(
		Prisma6GPCount,
		DenseMatrix(3, Prisma6NodesCount)
	);

	std::vector<std::vector<double> > rst = Prisma6_rst();

	for (unsigned int i = 0; i < Prisma6GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

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
		m(2, 2) =          - s / 2.0;
		m(2, 3) = -r / 2.0 - s / 2.0 + 1.0 / 2.0;
		m(2, 4) =  r / 2.0;
		m(2, 5) =            s / 2.0;
	}

	return dN;
}

static std::vector<DenseMatrix> Prisma6_N() {
	std::vector<DenseMatrix> N(
		Prisma6GPCount,
		DenseMatrix(1, Prisma6NodesCount)
	);

	std::vector<std::vector<double> > rst = Prisma6_rst();

	for (unsigned int i = 0; i < Prisma6GPCount; i++) {
		DenseMatrix &m = N[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

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

static std::vector<double> Prisma6_weight()
{
	std::vector<double> w;
	switch (Prisma6GPCount) {
	case 9: {
		double v1 = 5.0 / 54.0;
		double v2 = 8.0 / 54.0;
		w = { v1, v1, v1, v2, v2, v2, v1, v1, v1 };
		break;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Prisma6 GP count.";
		exit(EXIT_FAILURE);
	}

	return w;
}

std::vector<DenseMatrix> Prisma6::_dN = Prisma6_dN();
std::vector<DenseMatrix> Prisma6::_N = Prisma6_N();
std::vector<double> Prisma6::_weighFactor = Prisma6_weight();

bool Prisma6::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Prisma6 is 3D element
	return false;
#endif

	switch (n) {
	case Prisma6NodesCount:
		for (eslocal i = 0; i < Prisma6NodesCount - 1; i++) {
			for (eslocal j = i + 1; j < Prisma6NodesCount; j++) {
				if (Element::match(indices, i, j)) {
					return false;
				}
			}
		}
		return true;
	case 8 : {
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
	default:
		return false;
	}
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

size_t Prisma6::fillEdges()
{
	eslocal line[Line2NodesCount];

	if (_edges.size() == Prisma6EdgeCount) {
		return Prisma6EdgeCount;
	}
	_edges.reserve(Prisma6EdgeCount);

	size_t filled = _edges.size();

	for (size_t e = 0 ; e < Prisma6EdgeCount; e++) {
		for (size_t n = 0; n < Line2NodesCount; n++) {
			line[n] = _indices[_edgesNodes[e][n]];
		}
		addUniqueEdge<Line2>(line, filled, Line2NodesCount);
	}

	return filled;
}

size_t Prisma6::fillFaces()
{
	eslocal square[Square4NodesCount];
	eslocal triangle[Triangle3NodesCount];
	if (_faces.size() == Prisma6FacesCount) {
		return Prisma6FacesCount;
	}
	_faces.reserve(Prisma6FacesCount);

	size_t filled = _faces.size();

	for (size_t f = 0; f < 3; f++) {
		for (size_t n = 0; n < Square4NodesCount; n++) {
			square[n] = _indices[_facesNodes[f][n]];
		}
		addUniqueFace<Square4>(square, filled, Square4NodesCount);
	}

	for (size_t f = 3 ; f < Prisma6FacesCount; f++) {
		for (size_t n = 0; n < Triangle3NodesCount; n++) {
			triangle[n] = _indices[_facesNodes[f][n]];
		}
		addUniqueFace<Triangle3>(triangle, filled, Triangle3NodesCount);
	}

	return filled;
}

Element* Prisma6::addFace(const std::vector<eslocal> &nodes)
{
	for (size_t f = 0; f < faces(); f++) {
		size_t found;
		for (found = 0; found < _facesNodes[f].size(); found++) {
			if (!std::binary_search(nodes.begin(), nodes.end(), _indices[_facesNodes[f][found]])) {
				break;
			}
		}
		if (found == _facesNodes[f].size()) {
			if (f < 3) {
				eslocal square[Square4NodesCount];
				for (size_t n = 0; n < Square4NodesCount; n++) {
					square[n] = _indices[_facesNodes[f][n]];
				}
				return addUniqueFace<Square4>(square, _faces.size(), Square4NodesCount);
			} else {
				eslocal triangle[Triangle3NodesCount];
				for (size_t n = 0; n < Triangle3NodesCount; n++) {
					triangle[n] = _indices[_facesNodes[f][n]];
				}
				return addUniqueFace<Triangle3>(triangle, _faces.size(), Triangle3NodesCount);
			}
		}
	}
	return NULL;
}

Prisma6::Prisma6(const eslocal *indices, eslocal n, const eslocal *params)
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
		ESINFO(ERROR) << "It is not possible to create Prisma6 from " << n << " indices.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Prisma6::Prisma6(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}






