#include "square4.h"

using namespace espreso;

size_t Square4::_counter = 0;

std::vector<Property> Square4::_DOFElement;
std::vector<Property> Square4::_DOFFace;
std::vector<Property> Square4::_DOFEdge;
std::vector<Property> Square4::_DOFPoint;
std::vector<Property> Square4::_DOFMidPoint;

static std::vector<DenseMatrix> get_dN() {
	std::vector<DenseMatrix> dN(
		Square4GPCount,
		DenseMatrix(2, Square4NodesCount)
	);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Square4GPCount; i++) {
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = 0.25 * ( t - 1);
		m(0, 1) = 0.25 * (-t + 1);
		m(0, 2) = 0.25 * ( t + 1);
		m(0, 3) = 0.25 * (-t - 1);

		// dNt - derivation of basis function
		m(1, 0) = 0.25 * ( s - 1);
		m(1, 1) = 0.25 * (-s - 1);
		m(1, 2) = 0.25 * ( s + 1);
		m(1, 3) = 0.25 * (-s + 1);

//		m(2, 0) = 0;
//		m(2, 1) = 0;
//		m(2, 2) = 0;
//		m(2, 3) = 0;
	}

	return dN;
}

static std::vector<DenseMatrix> get_N() {
	std::vector<DenseMatrix> N(
		Square4GPCount,
		DenseMatrix(1, Square4NodesCount)
	);

	double CsQ_scale = 0.5773502691896258;

	for (unsigned int i = 0; i < Square4GPCount; i++) {
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		N[i](0, 0) = 0.25 * (1 - s) * (1 - t);
		N[i](0, 1) = 0.25 * (s + 1) * (1 - t);
		N[i](0, 2) = 0.25 * (s + 1) * (t + 1);
		N[i](0, 3) = 0.25 * (1 - s) * (t + 1);
	}

	return N;
}

std::vector<DenseMatrix> Square4::_dN = get_dN();
std::vector<DenseMatrix> Square4::_N = get_N();
std::vector<double> Square4::_weighFactor(Square4GPCount, 1);

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

std::vector<eslocal> Square4::getFace(size_t face) const
{
	if (face < 3) {
		return std::vector<eslocal>(_indices + face, _indices + face + 2);
	} else {
		return std::vector<eslocal>({_indices[3], _indices[0]});
	}
}

Element* getF(const eslocal *indices, const eslocal *params, size_t face)
{
	std::vector<eslocal> result(2);
	if (face < 3) {
		result[0] = indices[face];
		result[1] = indices[face + 1];
	} else {
		result[0] = indices[3];
		result[1] = indices[0];
	}

	return new Line2(indices, params);
}

Element* Square4::getFullFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Element* Square4::getCoarseFace(size_t face) const
{
	ESINFO(ERROR) << "get FACE is not implemented";
	return NULL;
}

Square4::Square4(const eslocal *indices, const eslocal *params): Element(params)
{
	memcpy(_indices, indices, Square4NodesCount * sizeof(eslocal));
}



