
#include "line3.h"
#include "../../../basis/matrices/denseMatrix.h"

using namespace espreso;

std::vector<Property> Line3::_DOFElement;
std::vector<Property> Line3::_DOFFace;
std::vector<Property> Line3::_DOFEdge;
std::vector<Property> Line3::_DOFPoint;
std::vector<Property> Line3::_DOFMidPoint;

static std::vector<DenseMatrix> get_dN() {
	std::vector<DenseMatrix> dN(
			Line3GPCount,
		DenseMatrix(1, Line3NodesCount)
	);

	std::vector<double> s = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

	for (unsigned int i = 0; i < Line3GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = s[i] - 1 / 2.0;
		m(0, 1) = -2 * s[i];
		m(0, 2) = s[i] + 1 / 2.0;
	}

	return dN;
}

static std::vector<DenseMatrix> get_N() {
	std::vector<DenseMatrix> N(
		Line3GPCount,
		DenseMatrix(1, Line3NodesCount)
	);

	std::vector<double> s = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

	for (unsigned int i = 0; i < Line3GPCount; i++) {
		N[i](0, 0) = (1 / 2.0) * (s[i] - 1) * s[i];
		N[i](0, 1) = 1 - s[i] * s[i];
		N[i](0, 2) = (1 / 2.0) * (s[i] + 1) * s[i];
	}

	return N;
}

std::vector<DenseMatrix> Line3::_dN = get_dN();
std::vector<DenseMatrix> Line3::_N = get_N();
std::vector<double> Line3::_weighFactor = { 5 / 9.0, 8 / 9.0, 5 / 9.0 };

bool Line3::match(const eslocal *indices, eslocal n)
{
	if (n != 3) {
		return false;
	}

	if (Element::match(indices, 0, 1) || Element::match(indices, 0, 2) || Element::match(indices, 1, 2)) {
		return false;
	}

	return true;
}

std::vector<eslocal> Line3::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(1);

	result[0] = (nodeIndex == 0) ? _indices[1] : _indices[0];

	return result;
}

Line3::Line3(const eslocal *indices)
{
	memcpy(_indices, indices, Line3NodesCount * sizeof(eslocal));
}

Line3::Line3(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	// line has not params
}



