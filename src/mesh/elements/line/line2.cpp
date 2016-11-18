
#include "line2.h"

using namespace espreso;

std::vector<Property> Line2::_DOFElement;
std::vector<Property> Line2::_DOFFace;
std::vector<Property> Line2::_DOFEdge;
std::vector<Property> Line2::_DOFPoint;
std::vector<Property> Line2::_DOFMidPoint;

static std::vector<DenseMatrix> get_dN() {
	std::vector<DenseMatrix> dN(
			Line2GPCount,
		DenseMatrix(1, Line2NodesCount)
	);

	for (unsigned int i = 0; i < Line2GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = -1 / 2.0;
		m(0, 1) =  1 / 2.0;
	}

	return dN;
}

static std::vector<DenseMatrix> get_N() {
	std::vector<DenseMatrix> N(
		Line2GPCount,
		DenseMatrix(1, Line2NodesCount)
	);

	std::vector<double> s = { 1 / sqrt(3), -1 / sqrt(3) };

	for (unsigned int i = 0; i < Line2GPCount; i++) {
		N[i](0, 0) = (1 - s[i]) / 2.0;
		N[i](0, 1) = (1 + s[i]) / 2.0;
	}

	return N;
}

std::vector<DenseMatrix> Line2::_dN = get_dN();
std::vector<DenseMatrix> Line2::_N = get_N();
std::vector<double> Line2::_weighFactor = { 1, 1 };

bool Line2::match(const eslocal *indices, eslocal n)
{
	if (n != 2) {
		return false;
	}

	if (Element::match(indices, 0, 1)) {
		return false;
	}

	return true;
}

std::vector<eslocal> Line2::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(1);

	result[0] = (nodeIndex == 0) ? _indices[1] : _indices[0];

	return result;
}

Line2::Line2(const eslocal *indices)
{
	memcpy(_indices, indices, Line2NodesCount * sizeof(eslocal));
}

Line2::Line2(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	// line has not params
}


