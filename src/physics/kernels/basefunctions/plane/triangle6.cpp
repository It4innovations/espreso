
#include "triangle.h"
#include "triangle6.h"


#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Triangle6::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, w;
	if (!Triangle::gpw(order, r, s, w)) {
		eslog::internalFailure("cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(r.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(r.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		(*N)[i](0, 0) = (1.0 - r[i] - s[i]) * (1.0 - 2.0 * (r[i] + s[i]));
		(*N)[i](0, 1) = -(r[i]) * (1.0 - 2.0 * r[i]);
		(*N)[i](0, 2) = -(s[i]) * (1.0 - 2.0 * s[i]);
		(*N)[i](0, 3) = 4.0 * (r[i]) * (1.0 - r[i] - s[i]);
		(*N)[i](0, 4) = 4.0 * (r[i]) * (s[i]);
		(*N)[i](0, 5) = 4.0 * (s[i]) * (1.0 - r[i] - s[i]);
	}

	for (size_t i = 0; i < r.size(); i++) {
		///dN contains [dNs, dNt]
		MatrixDense &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -3.0 + 4.0 * r[i] + 4.0 * s[i];
		m(0, 1) = -1.0 + 4.0 * r[i];
		m(0, 2) = 0.0;
		m(0, 3) = 4.0 - 8.0 * r[i] - 4.0 * s[i];
		m(0, 4) = 4.0 * s[i];
		m(0, 5) = -4.0 * s[i];

		// dNt - derivation of basis function
		m(1, 0) = -3.0 + 4.0 * r[i] + 4.0 * s[i];
		m(1, 1) = 0.0;
		m(1, 2) = -1.0 + 4.0 * s[i];
		m(1, 3) = -4.0 * r[i];
		m(1, 4) = 4.0 * r[i];
		m(1, 5) = 4.0 - 4.0 * r[i] - 8.0 * s[i];
	}
}

void Triangle6::setBaseFunctions(Element &self)
{
	size_t GPCount = 6, nodeCount = 6;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(2, nodeCount));
	self.weighFactor = new std::vector<double>(
		{ 0.109951743655322 / 2.0, 0.109951743655322 / 2.0, 0.109951743655322 / 2.0,
		  0.223381589678011 / 2.0, 0.223381589678011 / 2.0, 0.223381589678011 / 2.0 });

	const std::vector<double> s = {
		0.091576213509771,
		0.816847572980459,
		0.091576213509771,
		0.445948490915965,
		0.108103018168070,
		0.445948490915965
	};
	const std::vector<double> t = {
		0.091576213509771,
		0.091576213509771,
		0.816847572980459,
		0.445948490915965,
		0.445948490915965,
		0.108103018168070
	};

	for (size_t i = 0; i < GPCount; i++) {
		(*self.N)[i](0, 0) = (1.0 - s[i] - t[i]) * (1.0 - 2.0 * (s[i] + t[i]));
		(*self.N)[i](0, 1) = -(s[i]) * (1.0 - 2.0 * s[i]);
		(*self.N)[i](0, 2) = -(t[i]) * (1.0 - 2.0 * t[i]);
		(*self.N)[i](0, 3) = 4.0 * (s[i]) * (1.0 - s[i] - t[i]);
		(*self.N)[i](0, 4) = 4.0 * (s[i]) * (t[i]);
		(*self.N)[i](0, 5) = 4.0 * (t[i]) * (1.0 - s[i] - t[i]);
	}

	for (size_t i = 0; i < GPCount; i++) {
		///dN contains [dNs, dNt]
		MatrixDense &m = (*self.dN)[i];

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

	BaseFunctions::created(self);
}
