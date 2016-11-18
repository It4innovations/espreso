
#include "triangle3.h"
#include "../line/line2.h"

using namespace espreso;

std::vector<Property> Triangle3::_DOFElement;
std::vector<Property> Triangle3::_DOFFace;
std::vector<Property> Triangle3::_DOFEdge;
std::vector<Property> Triangle3::_DOFPoint;
std::vector<Property> Triangle3::_DOFMidPoint;

static std::vector<std::vector<double> > get_st()
{
	std::vector< std::vector<double> > st(2, std::vector<double>(Triangle3GPCount));

	switch (Triangle3GPCount) {
	case 1:
		st[0] = { 1.0 / 3 };
		st[1] = { 1.0 / 3 };
		return st;
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
		Triangle3GPCount,
		DenseMatrix(2, Triangle3NodesCount)
	);

	std::vector<std::vector<double> > st = get_st();
	const std::vector<double> &s = st[0];
	const std::vector<double> &t = st[1];

	for (unsigned int i = 0; i < Triangle3GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = -1;
		m(0, 1) =  1;
		m(0, 2) =  0;

		// dNt - derivation of basis function
		m(1, 0) = -1;
		m(1, 1) =  0;
		m(1, 2) =  1;
	}

	return dN;
}

static std::vector<DenseMatrix> get_N()
{
	std::vector<DenseMatrix> N(
		Triangle3GPCount,
		DenseMatrix(1, Triangle3NodesCount)
	);

	std::vector<std::vector<double> > st = get_st();
	const std::vector<double> &s = st[0];
	const std::vector<double> &t = st[1];

	for (unsigned int i = 0; i < Triangle3GPCount; i++) {
		N[i](0, 0) = 1 - s[i] - t[i];
		N[i](0, 1) = s[i];
		N[i](0, 2) = t[i];
	}

	return N;
}

static std::vector<double> get_w()
{
	switch (Triangle3GPCount) {
	case 1:
		return { 1. / 2 };
	case 6:
		return { 0.109951743655322 / 2.0, 0.109951743655322 / 2.0,  0.109951743655322 / 2.0, 0.223381589678011 / 2.0, 0.223381589678011 / 2.0, 0.223381589678011 / 2.0 };
	default:
		ESINFO(ERROR) << "Unknown number of Triangle3 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Triangle3::_dN = get_dN();
std::vector<DenseMatrix> Triangle3::_N = get_N();
std::vector<double> Triangle3::_weighFactor = get_w();

bool Triangle3::match(eslocal *indices, eslocal n)
{
	if (n != 3) {
		return false;
	}

	for (eslocal i = 0; i < Triangle3NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Triangle3NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Triangle3::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	result.reserve(2);

	for (eslocal i = 0; i < Triangle3NodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}

	return result;
}

void Triangle3::fillEdges()
{
	eslocal line[Line2NodesCount];
	_edges.reserve(Triangle3EdgeCount);

	size_t filled = _edges.size();

	for (size_t edge = 0; edge < 3; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 3];
		addEdge<Line2>(_edges, line, filled);
	}
}

Triangle3::Triangle3(const eslocal *indices)
{
	memcpy(_indices, indices, Triangle3NodesCount * sizeof(eslocal));
}

Triangle3::Triangle3(const eslocal *indices, const eslocal *params)
{
	memcpy(_indices, indices, Triangle3NodesCount * sizeof(eslocal));
	_params.insert(_params.end(), params, params + PARAMS_SIZE);
}

Triangle3::Triangle3(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	if (params) {
		_params.resize(params);
		is.read(reinterpret_cast<char *>(_params.data()), sizeof(eslocal) * params);
	}
}

