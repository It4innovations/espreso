
#include "tetrahedron.h"
#include "tetrahedron10.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Tetrahedron10::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, t, w;
	if (!Tetrahedron::gpw(order, r, s, t, w)) {
		eslog::globalerror("ESPRESO internal error: cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(s.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(s.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*N)[i];

		m(0, 0) = r[i] * (2.0 * r[i] - 1.0);
		m(0, 1) = s[i] * (2.0 * s[i] - 1.0);
		m(0, 2) = t[i] * (2.0 * t[i] - 1.0);
		m(0, 3) = 2.0 * r[i] * r[i] + 4.0 * r[i] * s[i] + 4.0 * r[i] * t[i] - 3.0 * r[i] + 2.0* s[i] * s[i] + 4.0 * s[i] * t[i] - 3.0 * s[i] + 2.0 * t[i] * t[i] - 3.0 * t[i] + 1.0;
		m(0, 4) = 4.0 * r[i] * s[i];
		m(0, 5) = 4.0 * s[i] * t[i];
		m(0, 6) = 4.0 * r[i] * t[i];
		m(0, 7) = r[i] * (-4.0 * r[i] - 4.0 * s[i] - 4.0 * t[i] + 4.0);
		m(0, 8) = s[i] * (-4.0 * r[i] - 4.0 * s[i] - 4.0 * t[i] + 4.0);
		m(0, 9) = t[i] * (-4.0 * r[i] - 4.0 * s[i] - 4.0 * t[i] + 4.0);
	}

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*dN)[i];

		m(0, 0) = 4.0 * r[i] - 1.0;
		m(0, 1) = 0;
		m(0, 2) = 0;
		m(0, 3) = 4.0 * r[i] + 4.0 * s[i] + 4.0 * t[i] - 3.0;
		m(0, 4) = 4.0 * s[i];
		m(0, 5) = 0;
		m(0, 6) = 4.0 * t[i];
		m(0, 7) = -8.0 * r[i] - 4.0 * s[i] - 4.0 * t[i] + 4.0;
		m(0, 8) = -4.0 * s[i];
		m(0, 9) = -4.0 * t[i];

		m(2, 0) = 0;
		m(2, 1) = 4.0 * s[i] - 1.0;
		m(2, 2) = 0 ;
		m(2, 3) = 4.0 * r[i] + 4.0 * s[i] + 4.0 * t[i] - 3.0;
		m(2, 4) = 4.0 * r[i];
		m(2, 5) = 4.0 * t[i];
		m(2, 6) = 0;
		m(2, 7) = -4.0 * r[i];
		m(2, 8) = -4.0 * r[i] - 8.0 * s[i] - 4.0 * t[i] + 4.0;
		m(2, 9) = -4.0 * t[i];

		m(1, 0) = 0;
		m(1, 1) = 0;
		m(1, 2) = 4.0 * t[i] - 1.0;
		m(1, 3) = 4.0 * r[i] + 4.0 * s[i] + 4.0* t[i]  - 3.0;
		m(1, 4) = 0;
		m(1, 5) = 4.0 * s[i];
		m(1, 6) = 4.0 * r[i];
		m(1, 7) = -4.0 * r[i];
		m(1, 8) = -4.0 * s[i];
		m(1, 9) = -4.0 * r[i] - 4.0 * s[i] - 8.0 * t[i] + 4.0;
	}
}

void Tetrahedron10::setBaseFunctions(Element &self)
{
	size_t GPCount = 15, nodeCount = 10;

	std::vector<std::vector<double> > rst(3);

	switch (GPCount) {
	case 15: {
		rst[0] = {
				0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
				0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
				0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
				0.0665501535736643, 0.4334498464263357, 0.4334498464263357};
		rst[1] = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
				0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
				0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
				0.4334498464263357, 0.0665501535736643, 0.4334498464263357};
		rst[2] = {
				0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
				0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
				0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
				0.4334498464263357, 0.4334498464263357, 0.0665501535736643};
		break;
	}
	default:
		exit(1);
	}

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(3, nodeCount));
	self.weighFactor = new std::vector<double>();

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		MatrixDense &m = (*self.N)[i];

		m(0, 0) = r * (2.0 * r - 1.0);
		m(0, 1) = s * (2.0 * s - 1.0);
		m(0, 2) = t * (2.0 * t - 1.0);
		m(0, 3) = 2.0 * r * r + 4.0 * r * s + 4.0 * r * t - 3.0 * r + 2.0* s * s + 4.0 * s * t - 3.0 * s + 2.0 * t * t - 3.0 * t + 1.0;
		m(0, 4) = 4.0 * r * s;
		m(0, 5) = 4.0 * s * t;
		m(0, 6) = 4.0 * r * t;
		m(0, 7) = r * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
		m(0, 8) = s * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
		m(0, 9) = t * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		MatrixDense &m = (*self.dN)[i];

		m(0, 0) = 4.0 * r - 1.0;
		m(0, 1) = 0;
		m(0, 2) = 0;
		m(0, 3) = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
		m(0, 4) = 4.0 * s;
		m(0, 5) = 0;
		m(0, 6) = 4.0 * t;
		m(0, 7) = -8.0 * r - 4.0 * s - 4.0 * t + 4.0;
		m(0, 8) = -4.0 * s;
		m(0, 9) = -4.0 * t;

		m(2, 0) = 0;
		m(2, 1) = 4.0 * s - 1.0;
		m(2, 2) = 0 ;
		m(2, 3) = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
		m(2, 4) = 4.0 * r;
		m(2, 5) = 4.0 * t;
		m(2, 6) = 0;
		m(2, 7) = -4.0 * r;
		m(2, 8) = -4.0 * r - 8.0 * s - 4.0 * t + 4.0;
		m(2, 9) = -4.0 * t;

		m(1, 0) = 0;
		m(1, 1) = 0;
		m(1, 2) = 4.0 * t - 1.0;
		m(1, 3) = 4.0 * r + 4.0 * s + 4.0* t  - 3.0;
		m(1, 4) = 0;
		m(1, 5) = 4.0 * s;
		m(1, 6) = 4.0 * r;
		m(1, 7) = -4.0 * r;
		m(1, 8) = -4.0 * s;
		m(1, 9) = -4.0 * r - 4.0 * s - 8.0 * t + 4.0;
	}

	switch (GPCount) {
	case 15: {
		self.weighFactor->resize( 1, 0.030283678097089);
		self.weighFactor->resize( 5, 0.006026785714286);
		self.weighFactor->resize( 9, 0.011645249086029);
		self.weighFactor->resize(15, 0.010949141561386);
		break;
	}
	default:
		exit(1);
	}

	BaseFunctions::created(self);
}






