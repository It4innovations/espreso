
#include "square4.h"
#include "../line/line2.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

std::vector<Property> Square4::_DOFElement;
std::vector<Property> Square4::_DOFFace;
std::vector<Property> Square4::_DOFEdge;
std::vector<Property> Square4::_DOFPoint;
std::vector<Property> Square4::_DOFMidPoint;

static std::vector<DenseMatrix> get_dN(double CsQ_scale) {
	std::vector<DenseMatrix> dN(
		Square4GPCount,
		DenseMatrix(2, Square4NodesCount)
	);

	std::vector<double> s = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
	std::vector<double> t = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };
	for (unsigned int i = 0; i < Square4GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = 0.25 * ( t[i] - 1);
		m(0, 1) = 0.25 * (-t[i] + 1);
		m(0, 2) = 0.25 * ( t[i] + 1);
		m(0, 3) = 0.25 * (-t[i] - 1);

		// dNt - derivation of basis function
		m(1, 0) = 0.25 * ( s[i] - 1);
		m(1, 1) = 0.25 * (-s[i] - 1);
		m(1, 2) = 0.25 * ( s[i] + 1);
		m(1, 3) = 0.25 * (-s[i] + 1);

//		m(2, 0) = 0;
//		m(2, 1) = 0;
//		m(2, 2) = 0;
//		m(2, 3) = 0;
	}

	return dN;
}

static std::vector<DenseMatrix> get_N(double CsQ_scale) {
	std::vector<DenseMatrix> N(
		Square4GPCount,
		DenseMatrix(1, Square4NodesCount)
	);

	std::vector<double> s = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
	std::vector<double> t = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };
	for (unsigned int i = 0; i < Square4GPCount; i++) {
		// basis function
		N[i](0, 0) = 0.25 * (1 - s[i]) * (1 - t[i]);
		N[i](0, 1) = 0.25 * (s[i] + 1) * (1 - t[i]);
		N[i](0, 2) = 0.25 * (s[i] + 1) * (t[i] + 1);
		N[i](0, 3) = 0.25 * (1 - s[i]) * (t[i] + 1);
	}

	return N;
}

std::vector<DenseMatrix> Square4::_dN = get_dN(0.577350269189626);
std::vector<DenseMatrix> Square4::_N = get_N(0.577350269189626);
std::vector<double> Square4::_weighFactor(Square4GPCount, 1);

std::vector<DenseMatrix> Square4::_v_dN = get_dN(1);
std::vector<DenseMatrix> Square4::_v_N = get_N(1);
std::vector<double> Square4::_v_weighFactor(Square4GPCount, 1);

bool Square4::match(const eslocal *indices, eslocal n)
{
	if (n != 4) {
		return false;
	}

	for (eslocal i = 0; i < Square4NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Square4NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Square4::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(2);

	result[0] = _indices[(nodeIndex + 1) % 4];
	result[1] = _indices[(nodeIndex + 3) % 4];

	return result;
}

size_t Square4::fillEdges()
{
	eslocal line[Line2NodesCount];

	if (_edges.size() == Square4EdgeCount) {
		return Square4EdgeCount;
	}
	_edges.reserve(Square4EdgeCount);

	size_t filled = _edges.size();

	for (size_t edge = 0; edge < 4; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 4];
		addUniqueEdge<Line2>(line, filled);
	}
	return filled;
}

Square4::Square4(const eslocal *indices)
{
	memcpy(_indices, indices, Square4NodesCount * sizeof(eslocal));
}

Square4::Square4(const eslocal *indices, const eslocal *params)
{
	memcpy(_indices, indices, Square4NodesCount * sizeof(eslocal));
	_params.insert(_params.end(), params, params + PARAMS_SIZE);
}

Square4::Square4(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	if (params) {
		_params.resize(params);
		is.read(reinterpret_cast<char *>(_params.data()), sizeof(eslocal) * params);
	}
}

