
#include "pyramid.h"
#include "pyramid5.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Pyramid5::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, t, w;
	if (!Pyramid::gpw(order, r, s, t, w)) {
		eslog::internalFailure("cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(s.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(s.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*N)[i];

		m(0, 0) = 0.125 * ((1 - r[i]) * (1 - s[i]) * (1 - t[i]));
		m(0, 1) = 0.125 * ((1 + r[i]) * (1 - s[i]) * (1 - t[i]));
		m(0, 2) = 0.125 * ((1 + r[i]) * (1 + s[i]) * (1 - t[i]));
		m(0, 3) = 0.125 * ((1 - r[i]) * (1 + s[i]) * (1 - t[i]));
		m(0, 4) = 0.125 * ( 4 * (1 + t[i]));
	}

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*dN)[i];

		m(0, 0) = 0.125 * (-(1. - s[i]) * (1. - t[i]));
		m(0, 1) = 0.125 * ( (1. - s[i]) * (1. - t[i]));
		m(0, 2) = 0.125 * ( (1. + s[i]) * (1. - t[i]));
		m(0, 3) = 0.125 * (-(1. + s[i]) * (1. - t[i]));
		m(0, 4) = 0;

		m(1, 0) = 0.125 * (-(1. - r[i]) * (1. - t[i]));
		m(1, 1) = 0.125 * (-(1. + r[i]) * (1. - t[i]));
		m(1, 2) = 0.125 * ( (1. + r[i]) * (1. - t[i]));
		m(1, 3) = 0.125 * ( (1. - r[i]) * (1. - t[i]));
		m(1, 4) = 0;

		m(2, 0) = 0.125 * (-(1. - r[i]) * (1. - s[i]));
		m(2, 1) = 0.125 * (-(1. + r[i]) * (1. - s[i]));
		m(2, 2) = 0.125 * (-(1. + r[i]) * (1. + s[i]));
		m(2, 3) = 0.125 * (-(1. - r[i]) * (1. + s[i]));
		m(2, 4) = 0.125 * (4.0);
	}
}

void Pyramid5::setBaseFunctions(Element &self)
{
	size_t GPCount = 8, nodeCount = 5;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(3, nodeCount));
	self.weighFactor = new std::vector<double>();

	std::vector< std::vector<double> > rst(3, std::vector<double>(GPCount));

	switch (GPCount) {
		case 8: {
			double v = 0.577350269189625953;
			rst[0] = {  v,  v,  v,  v, -v, -v, -v, -v };
			rst[1] = { -v, -v,  v,  v, -v, -v,  v,  v };
			rst[2] = { -v,  v, -v,  v, -v,  v, -v,  v };
			break;
		}
		default:
			exit(1);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		MatrixDense &m = (*self.N)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// basis function
		m(0, 0) = 0.125 * ((1 - r) * (1 - s) * (1 - t));
		m(0, 1) = 0.125 * ((1 + r) * (1 - s) * (1 - t));
		m(0, 2) = 0.125 * ((1 + r) * (1 + s) * (1 - t));
		m(0, 3) = 0.125 * ((1 - r) * (1 + s) * (1 - t));
		m(0, 4) = 0.125 * ( 4 * (1 + t));
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		MatrixDense &m = (*self.dN)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// dNr - derivation of basis function
		m(0, 0) = 0.125 * (-(1. - s) * (1. - t));
		m(0, 1) = 0.125 * ( (1. - s) * (1. - t));
		m(0, 2) = 0.125 * ( (1. + s) * (1. - t));
		m(0, 3) = 0.125 * (-(1. + s) * (1. - t));
		m(0, 4) = 0;

		// dNs - derivation of basis function
		m(1, 0) = 0.125 * (-(1. - r) * (1. - t));
		m(1, 1) = 0.125 * (-(1. + r) * (1. - t));
		m(1, 2) = 0.125 * ( (1. + r) * (1. - t));
		m(1, 3) = 0.125 * ( (1. - r) * (1. - t));
		m(1, 4) = 0;

		// dNt - derivation of basis function
		m(2, 0) = 0.125 * (-(1. - r) * (1. - s));
		m(2, 1) = 0.125 * (-(1. + r) * (1. - s));
		m(2, 2) = 0.125 * (-(1. + r) * (1. + s));
		m(2, 3) = 0.125 * (-(1. - r) * (1. + s));
		m(2, 4) = 0.125 * (4.0);
	}

	switch (GPCount) {
	case 8: {
		(*self.weighFactor) = std::vector<double> (8, 1.0);
		break;
	}
	default:
		exit(1);
	}

	BaseFunctions::created(self);
}

