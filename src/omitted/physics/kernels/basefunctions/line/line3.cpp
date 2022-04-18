
#include "line.h"
#include "line3.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Line3::setGaussPointsForOrder(int order)
{
	std::vector<double> r, w;
	if (!Line::gpw(order, r, w)) {
		eslog::internalFailure("cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(r.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(r.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		(*N)[i](0, 0) = (1 / 2.0) * (r[i] - 1) * r[i];
		(*N)[i](0, 1) = 1 - r[i] * r[i];
		(*N)[i](0, 2) = (1 / 2.0) * (r[i] + 1) * r[i];
	}

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = r[i] - 1 / 2.0;
		m(0, 1) = -2 * r[i];
		m(0, 2) = r[i] + 1 / 2.0;
	}

	BaseFunctions::created(*this);
}

void Line3::setBaseFunctions(Element &self)
{
	size_t GPCount = 3, nodeCount = 3;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.weighFactor = new std::vector<double>({ 5 / 9.0, 8 / 9.0, 5 / 9.0 });

	std::vector<double> s = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

	for (size_t i = 0; i < GPCount; i++) {
		(*self.N)[i](0, 0) = (1 / 2.0) * (s[i] - 1) * s[i];
		(*self.N)[i](0, 1) = (1 / 2.0) * (s[i] + 1) * s[i];
		(*self.N)[i](0, 2) = 1 - s[i] * s[i];
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		MatrixDense &m = (*self.dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = s[i] - 1 / 2.0;
		m(0, 1) = s[i] + 1 / 2.0;
		m(0, 2) = -2 * s[i];
	}

	BaseFunctions::created(self);
}



