
#include "square8.h"
#include "../line/line3.h"

using namespace espreso;

std::vector<Property> Square8::_DOFElement;
std::vector<Property> Square8::_DOFFace;
std::vector<Property> Square8::_DOFEdge;
std::vector<Property> Square8::_DOFPoint;
std::vector<Property> Square8::_DOFMidPoint;

static std::vector<std::vector<double> > get_st()
{
	std::vector< std::vector<double> > st(2, std::vector<double>(Square8GPCount));

	switch (Square8GPCount) {
	case 9: {
		double v = sqrt(0.6);
		st[0] = { -v,  v,  v, -v,  0,  v,  0, -v, 0 };
		st[1] = { -v, -v,  v,  v, -v,  0,  v,  0, 0 };
		return st;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Square8 GP count.";
		exit(EXIT_FAILURE);
	}
}

static std::vector<DenseMatrix> get_dN()
{
	std::vector<DenseMatrix> dN(
		Square8GPCount,
		DenseMatrix(2, Square8NodesCount)
	);

	std::vector<std::vector<double> > st = get_st();
	const std::vector<double> &s = st[0];
	const std::vector<double> &t = st[1];

	for (unsigned int i = 0; i < Square8GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = -((2 * s[i] + t[i]) * (t[i] - 1)) / 4;
		m(0, 1) = -((2 * s[i] - t[i]) * (t[i] - 1)) / 4;
		m(0, 2) = ((2 * s[i] + t[i]) * (t[i] + 1)) / 4;
		m(0, 3) = ((2 * s[i] - t[i]) * (t[i] + 1)) / 4;
		m(0, 4) = s[i] * (t[i] - 1);
		m(0, 5) = 1. / 2 - t[i] * t[i] / 2;
		m(0, 6) = -s[i] * (t[i] + 1);
		m(0, 7) = t[i] * t[i] / 2 - 1. / 2;

		// dNt - derivation of basis function
		m(1, 0) = -((s[i] + 2 * t[i]) * (s[i] - 1)) / 4;
		m(1, 1) = -((s[i] - 2 * t[i]) * (s[i] + 1)) / 4;
		m(1, 2) = ((s[i] + 2 * t[i]) * (s[i] + 1)) / 4;
		m(1, 3) = ((s[i] - 2 * t[i]) * (s[i] - 1)) / 4;
		m(1, 4) = s[i] * s[i] / 2 - 1. / 2;
		m(1, 5) = -t[i] * (s[i] + 1);
		m(1, 6) = 1. / 2 - s[i] * s[i] / 2;
		m(1, 7) = t[i] * (s[i] - 1);
	}

	return dN;
}

static std::vector<DenseMatrix> get_N() {
	std::vector<DenseMatrix> N(
		Square8GPCount,
		DenseMatrix(1, Square8NodesCount)
	);


	std::vector<std::vector<double> > st = get_st();
	const std::vector<double> &s = st[0];
	const std::vector<double> &t = st[1];

	for (unsigned int i = 0; i < Square8GPCount; i++) {
		N[i](0, 0) = -.25 * (s[i] - 1) * (t[i] - 1) * (s[i] + t[i] + 1);
		N[i](0, 1) =  .25 * (t[i] - 1) * (-s[i] * s[i] + t[i] * s[i] + t[i] + 1);
		N[i](0, 2) =  .25 * (s[i] + 1) * (t[i] + 1) * (s[i] + t[i] - 1);
		N[i](0, 3) =  .25 * (s[i] - 1) * (s[i] - t[i] + 1) * (t[i] + 1);
		N[i](0, 4) =  .5  * (s[i] * s[i] - 1) * (t[i] - 1);
		N[i](0, 5) = -.5  * (s[i] + 1) * (t[i] * t[i] - 1);
		N[i](0, 6) = -.5  * (s[i] * s[i] - 1) * (t[i] + 1);
		N[i](0, 7) =  .5  * (s[i] - 1) * (t[i] * t[i] - 1);
	}

	return N;
}

std::vector<DenseMatrix> Square8::_dN = get_dN();
std::vector<DenseMatrix> Square8::_N = get_N();
std::vector<double> Square8::_weighFactor = {
		25 / 81.0, 25 / 81.0, 25 / 81.0, 25 / 81.0,
		40 / 81.0, 40 / 81.0, 40 / 81.0, 40 / 81.0,
		64 / 81.0 };

bool Square8::match(const eslocal *indices, eslocal n)
{
	if (n != 8) {
		return false;
	}

	for (eslocal i = 0; i < Square8NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Square8NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Square8::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(2);

	if (nodeIndex < 4) {
		result[0] = _indices[nodeIndex + 4];
		result[1] = _indices[(nodeIndex + 3) % 4 + 4];
	} else {
		result[0] = _indices[(nodeIndex + 5) % 4];
		result[1] = _indices[nodeIndex - 4];
	}

	return result;
}

size_t Square8::fillEdges()
{
	eslocal line[Line3NodesCount];

	if (_edges.size() == Square8EdgeCount) {
		return Square8EdgeCount;
	}
	_edges.reserve(Square8EdgeCount);

	size_t filled = _edges.size();

	for (size_t edge = 0; edge < 4; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 4];
		line[2] = _indices[ edge + 4     ];
		addUniqueEdge<Line3>(line, filled);
	}
	return filled;
}

Square8::Square8(const eslocal *indices)
{
	memcpy(_indices, indices, Square8NodesCount * sizeof(eslocal));
}

Square8::Square8(const eslocal *indices, const eslocal *params)
{
	memcpy(_indices, indices, Square8NodesCount * sizeof(eslocal));
	_params.insert(_params.end(), params, params + PARAMS_SIZE);
}

Square8::Square8(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	if (params) {
		_params.resize(params);
		is.read(reinterpret_cast<char *>(_params.data()), sizeof(eslocal) * params);
	}
}


