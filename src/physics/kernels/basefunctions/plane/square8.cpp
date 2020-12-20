
#include "square.h"
#include "square8.h"


#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Square8::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, w;
	if (!Square::gpw(order, r, s, w)) {
		eslog::internalFailure("cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(r.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(r.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		(*N)[i](0, 0) = -.25 * (r[i] - 1) * (s[i] - 1) * (r[i] + s[i] + 1);
		(*N)[i](0, 1) =  .25 * (s[i] - 1) * (-r[i] * r[i] + s[i] * r[i] + s[i] + 1);
		(*N)[i](0, 2) =  .25 * (r[i] + 1) * (s[i] + 1) * (r[i] + s[i] - 1);
		(*N)[i](0, 3) =  .25 * (r[i] - 1) * (r[i] - s[i] + 1) * (s[i] + 1);
		(*N)[i](0, 4) =  .5  * (r[i] * r[i] - 1) * (s[i] - 1);
		(*N)[i](0, 5) = -.5  * (r[i] + 1) * (s[i] * s[i] - 1);
		(*N)[i](0, 6) = -.5  * (r[i] * r[i] - 1) * (s[i] + 1);
		(*N)[i](0, 7) =  .5  * (r[i] - 1) * (s[i] * s[i] - 1);
	}

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -((2 * r[i] + s[i]) * (s[i] - 1)) / 4;
		m(0, 1) = -((2 * r[i] - s[i]) * (s[i] - 1)) / 4;
		m(0, 2) = ((2 * r[i] + s[i]) * (s[i] + 1)) / 4;
		m(0, 3) = ((2 * r[i] - s[i]) * (s[i] + 1)) / 4;
		m(0, 4) = r[i] * (s[i] - 1);
		m(0, 5) = 1. / 2 - s[i] * s[i] / 2;
		m(0, 6) = -r[i] * (s[i] + 1);
		m(0, 7) = s[i] * s[i] / 2 - 1. / 2;

		// dNt - derivation of basis function
		m(1, 0) = -((r[i] + 2 * s[i]) * (r[i] - 1)) / 4;
		m(1, 1) = -((r[i] - 2 * s[i]) * (r[i] + 1)) / 4;
		m(1, 2) = ((r[i] + 2 * s[i]) * (r[i] + 1)) / 4;
		m(1, 3) = ((r[i] - 2 * s[i]) * (r[i] - 1)) / 4;
		m(1, 4) = r[i] * r[i] / 2 - 1. / 2;
		m(1, 5) = -s[i] * (r[i] + 1);
		m(1, 6) = 1. / 2 - r[i] * r[i] / 2;
		m(1, 7) = s[i] * (r[i] - 1);
	}

	BaseFunctions::created(*this);
}

void Square8::setBaseFunctions(Element &self)
{
	size_t GPCount = 9, nodeCount = 8;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(2, nodeCount));
	self.weighFactor = new std::vector<double>({
			25 / 81.0, 25 / 81.0, 25 / 81.0, 25 / 81.0,
			40 / 81.0, 40 / 81.0, 40 / 81.0, 40 / 81.0,
			64 / 81.0 });

	std::vector< std::vector<double> > st(2, std::vector<double>(GPCount));
	double v = sqrt(0.6);
	st[0] = { -v,  v,  v, -v,  0,  v,  0, -v, 0 };
	st[1] = { -v, -v,  v,  v, -v,  0,  v,  0, 0 };

	for (size_t i = 0; i < GPCount; i++) {
		const std::vector<double> &s = st[0];
		const std::vector<double> &t = st[1];

		(*self.N)[i](0, 0) = -.25 * (s[i] - 1) * (t[i] - 1) * (s[i] + t[i] + 1);
		(*self.N)[i](0, 1) =  .25 * (t[i] - 1) * (-s[i] * s[i] + t[i] * s[i] + t[i] + 1);
		(*self.N)[i](0, 2) =  .25 * (s[i] + 1) * (t[i] + 1) * (s[i] + t[i] - 1);
		(*self.N)[i](0, 3) =  .25 * (s[i] - 1) * (s[i] - t[i] + 1) * (t[i] + 1);
		(*self.N)[i](0, 4) =  .5  * (s[i] * s[i] - 1) * (t[i] - 1);
		(*self.N)[i](0, 5) = -.5  * (s[i] + 1) * (t[i] * t[i] - 1);
		(*self.N)[i](0, 6) = -.5  * (s[i] * s[i] - 1) * (t[i] + 1);
		(*self.N)[i](0, 7) =  .5  * (s[i] - 1) * (t[i] * t[i] - 1);
	}

	for (size_t i = 0; i < GPCount; i++) {
		MatrixDense &m = (*self.dN)[i];
		const std::vector<double> &s = st[0];
		const std::vector<double> &t = st[1];

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

	BaseFunctions::created(self);
}



