
#include "triangle6.h"
#include "../line/line3.h"

using namespace espreso;

std::vector<Property> Triangle6::_DOFElement;
std::vector<Property> Triangle6::_DOFFace;
std::vector<Property> Triangle6::_DOFEdge;
std::vector<Property> Triangle6::_DOFPoint;
std::vector<Property> Triangle6::_DOFMidPoint;

static std::vector<std::vector<double> > get_st()
{
	std::vector< std::vector<double> > st(2, std::vector<double>(Triangle6GPCount));

	switch (Triangle6GPCount) {
	case 6: {
		st[0] = {
			0.091576213509771,
			0.816847572980459,
			0.091576213509771,
			0.445948490915965,
			0.108103018168070,
			0.445948490915965
		};

		st[1] =  {
			0.091576213509771,
			0.091576213509771,
			0.816847572980459,
			0.445948490915965,
			0.445948490915965,
			0.108103018168070
		};
		return st;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Triangle3 GP count.";
		exit(EXIT_FAILURE);
	}
}

static std::vector<DenseMatrix> get_dN()
{
	std::vector<DenseMatrix> dN(
			Triangle6GPCount,
		DenseMatrix(2, Triangle6NodesCount)
	);

	std::vector<std::vector<double> > st = get_st();
	const std::vector<double> &s = st[0];
	const std::vector<double> &t = st[1];

	for (unsigned int i = 0; i < Triangle6GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = -3.0 + 4.0 * s[i] + 4.0 * t[i];
		m(0, 1) = -1.0 + 4.0 * s[i];
		m(0, 2) = 0.0;
		m(0, 3) = 4.0 - 8.0 * s[i] - 4.0 * t[i];
		m(0, 4) = 4.0 * t[i];
		m(0, 5) = -4.0 * t[i];

		// dNt - derivation of basis function
		m(1, 0) = -3.0 + 4.0 * s[i] + 4.0 * t[i];
		m(1, 1) = 0.0;
		m(1, 2) = -1.0 + 4.0 * t[i];
		m(1, 3) = -4.0 * s[i];
		m(1, 4) = 4.0 * s[i];
		m(1, 5) = 4.0 - 4.0 * s[i] - 8.0 * t[i];
	}

	return dN;
}

static std::vector<DenseMatrix> get_N() {
	std::vector<DenseMatrix> N(
			Triangle6GPCount,
		DenseMatrix(1, Triangle6NodesCount)
	);

	std::vector<std::vector<double> > st = get_st();
	const std::vector<double> &s = st[0];
	const std::vector<double> &t = st[1];

	for (unsigned int i = 0; i < Triangle6GPCount; i++) {
		N[i](0, 0) = (1.0 - s[i] - t[i]) * (1.0 - 2.0 * (s[i] + t[i]));
		N[i](0, 1) = -(s[i]) * (1.0 - 2.0 * s[i]);
		N[i](0, 2) = -(t[i]) * (1.0 - 2.0 * t[i]);
		N[i](0, 3) = 4.0 * (s[i]) * (1.0 - s[i] - t[i]);
		N[i](0, 4) = 4.0 * (s[i]) * (t[i]);
		N[i](0, 5) = 4.0 * (t[i]) * (1.0 - s[i] - t[i]);
	}

	return N;
}

std::vector<DenseMatrix> Triangle6::_dN = get_dN();
std::vector<DenseMatrix> Triangle6::_N = get_N();
std::vector<double> Triangle6::_weighFactor = { 0.109951743655322 / 2.0, 0.109951743655322 / 2.0,  0.109951743655322 / 2.0, 0.223381589678011 / 2.0, 0.223381589678011 / 2.0, 0.223381589678011 / 2.0 };

bool Triangle6::match(const eslocal *indices, const eslocal n)
{
	if (n != 6) {
		return false;
	}

	for (eslocal i = 0; i < Triangle6NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Triangle6NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Triangle6::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(2);

	if (nodeIndex < 3) {
		result[0] = _indices[nodeIndex + 3];
		result[1] = _indices[(nodeIndex + 2) % 3 + 3];
	} else {
		result[0] = _indices[(nodeIndex + 1) % 3];
		result[1] = _indices[nodeIndex - 3];
	}

	return result;
}

size_t Triangle6::fillEdges()
{
	eslocal line[Line3NodesCount];

	if (_edges.size() == Triangle6EdgeCount) {
		return Triangle6EdgeCount;
	}
	_edges.reserve(Triangle6EdgeCount);

	size_t filled = _edges.size();

	for (size_t edge = 0; edge < 3; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 3];
		line[2] = _indices[ edge + 3     ];
		addUniqueEdge<Line3>(line, filled);
	}
	return filled;
}

Triangle6::Triangle6(const eslocal *indices)
{
	memcpy(_indices, indices, Triangle6NodesCount * sizeof(eslocal));
}

Triangle6::Triangle6(const eslocal *indices, const eslocal *params)
{
	memcpy(_indices, indices, Triangle6NodesCount * sizeof(eslocal));
	_params.insert(_params.end(), params, params + PARAMS_SIZE);
}

Triangle6::Triangle6(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	if (params) {
		_params.resize(params);
		is.read(reinterpret_cast<char *>(_params.data()), sizeof(eslocal) * params);
	}
}

