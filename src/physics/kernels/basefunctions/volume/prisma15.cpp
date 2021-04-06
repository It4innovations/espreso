
#include "prisma.h"
#include "prisma15.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Prisma15::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, t, w;
	if (!Prisma::gpw(order, r, s, t, w)) {
		eslog::globalerror("ESPRESO internal error: cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(s.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(s.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*N)[i];

		// bas[i]is[i] funct[i]ion
		m(0, 0) = -(1.0 - r[i] - s[i]) * (1.0 - t[i]) * (2.0 * r[i] + 2.0 * s[i] + t[i]) / 2.0;
		m(0, 1) = r[i] * (1.0 - t[i]) * (2.0 * r[i] - t[i] - 2.0) / 2.0;
		m(0, 2) = s[i] * (1.0 - t[i]) * (2.0 * s[i] - t[i] - 2.0) / 2.0;
		m(0, 3) = -(1.0 - r[i] - s[i]) * (1.0 + t[i]) * (2.0 * r[i] + 2.0 * s[i] - t[i]) / 2.0;
		m(0, 4) = r[i] * (t[i] + 1.0) * (2.0 * r[i] + t[i] - 2.0) / 2.0;
		m(0, 5) = s[i] * (t[i] + 1.0) * (2.0 * s[i] + t[i] - 2.0) / 2.0;
		m(0, 6) = 2.0 * r[i] * (1.0 - r[i] - s[i]) * (1.0 - t[i]);
		m(0, 7) = 2.0 * r[i] * s[i] * (1.0 - t[i]);
		m(0, 8) = 2.0 * s[i] * (1.0 - r[i] - s[i]) * (1.0 - t[i]);
		m(0, 9) = 2.0 * r[i] * (1.0 - r[i] - s[i]) * (1.0 + t[i]);
		m(0, 10) = 2.0 * r[i] * s[i] * (1.0 + t[i]);
		m(0, 11) = 2.0 * s[i] * (1.0 - r[i] - s[i]) * (1.0 + t[i]);
		m(0, 12) = (1.0 - r[i] - s[i]) * (1.0 - pow(t[i], 2.0));
		m(0, 13) = r[i] * (1.0 - pow(t[i], 2.0));
		m(0, 14) = s[i] * (1.0 - pow(t[i], 2.0));
	}

	for (size_t i = 0; i < r.size(); i++) {
		///dN cont[i]ains[i] [dNr[i], dNs[i], dNt[i]]
		MatrixDense &m = (*dN)[i];

		// dNr[i] - der[i]ivat[i]ion of bas[i]is[i] funct[i]ion
		m(0, 0) = -(t[i] - 1.0) * (r[i] + s[i] - 1.0) - ((t[i] - 1.0) * (2.0 * r[i] + 2.0 * s[i] + t[i])) / 2.0;
		m(0, 1) = ((t[i] - 1.0) * (t[i] - 2.0 * r[i] + 2.0)) / 2.0 - r[i] * (t[i] - 1.0);
		m(0, 2) = 0.0;
		m(0, 3) = ((t[i] + 1.0) * (2.0 * r[i] + 2.0 * s[i] - t[i])) / 2.0 + (t[i] + 1.0) * (r[i] + s[i] - 1.0);
		m(0, 4) = r[i] * (t[i] + 1.0) + ((t[i] + 1.0) * (2.0 * r[i] + t[i] - 2.0)) / 2.0;
		m(0, 5) = 0.0;
		m(0, 6) = 2.0 * (t[i] - 1.0) * (r[i] + s[i] - 1.0) + 2.0 * r[i] * (t[i] - 1.0);
		m(0, 7) = (-2.0) * s[i] * (t[i] - 1.0);
		m(0, 8) = 2.0 * s[i] * (t[i] - 1.0);
		m(0, 9) = -2.0 * (t[i] + 1.0) * (r[i] + s[i] - 1.0) - 2.0 * r[i] * (t[i] + 1.0);
		m(0, 10) = 2.0 * s[i] * (t[i] + 1.0);
		m(0, 11) = -2.0 * s[i] * (t[i] + 1.0);
		m(0, 12) = pow(t[i], 2.0) - 1.0;
		m(0, 13) = 1.0 - pow(t[i], 2.0);
		m(0, 14) = 0.0;

		// dNs[i] - der[i]ivat[i]ion of bas[i]is[i] funct[i]ion
		m(1, 0) = -(t[i] - 1.0) * (r[i] + s[i] - 1.0) - ((t[i] - 1.0) * (2.0 * r[i] + 2.0 * s[i] + t[i])) / 2.0;
		m(1, 1) = 0.0;
		m(1, 2) = ((t[i] - 1.0) * (t[i] - 2.0 * s[i] + 2.0)) / 2.0 - s[i] * (t[i] - 1.0);
		m(1, 3) = ((t[i] + 1.0) * (2.0 * r[i] + 2.0 * s[i] - t[i])) / 2.0 + (t[i] + 1.0) * (r[i] + s[i] - 1.0);
		m(1, 4) = 0.0;
		m(1, 5) = s[i] * (t[i] + 1.0) + ((t[i] + 1.0) * (2.0 * s[i] + t[i] - 2.0)) / 2.0;
		m(1, 6) = 2.0 * r[i] * (t[i] - 1.0);
		m(1, 7) = (-2.0) * r[i] * (t[i] - 1.0);
		m(1, 8) = 2.0 * (t[i] - 1.0) * (r[i] + s[i] - 1.0) + 2.0 * s[i] * (t[i] - 1.0);
		m(1, 9) = (-2.0) * r[i] * (t[i] + 1.0);
		m(1, 10) = 2.0 * r[i] * (t[i] + 1.0);
		m(1, 11) = -2.0 * (t[i] + 1.0) * (r[i] + s[i] - 1.0) - 2.0 * s[i] * (t[i] + 1.0);
		m(1, 12) = pow(t[i], 2.0) - 1.0;
		m(1, 13) = 0.0;
		m(1, 14) = 1.0 - pow(t[i], 2.0);

		// dNt[i] - der[i]ivat[i]ion of bas[i]is[i] funct[i]ion
		m(2, 0) = -((r[i] + s[i] - 1.0) * (2.0 * r[i] + 2.0 * s[i] + t[i])) / 2.0 - ((t[i] - 1.0) * (r[i] + s[i] - 1.0)) / 2.0;
		m(2, 1) = (r[i] * (t[i] - 2.0 * r[i] + 2.0)) / 2.0 + (r[i] * (t[i] - 1.0)) / 2.0;
		m(2, 2) = (s[i] * (t[i] - 2.0 * s[i] + 2.0)) / 2.0 + (s[i] * (t[i] - 1.0)) / 2.0;
		m(2, 3) = ((r[i] + s[i] - 1.0) * (2.0 * r[i] + 2.0 * s[i] - t[i])) / 2.0 - ((t[i] + 1.0) * (r[i] + s[i] - 1.0)) / 2.0;
		m(2, 4) = (r[i] * (2.0 * r[i] + t[i] - 2.0)) / 2.0 + (r[i] * (t[i] + 1.0)) / 2.0;
		m(2, 5) = (s[i] * (2.0 * s[i] + t[i] - 2.0)) / 2.0 + (s[i] * (t[i] + 1.0)) / 2.0;
		m(2, 6) = 2.0 * r[i] * (r[i] + s[i] - 1.0);
		m(2, 7) = (-2.0) * r[i] * s[i];
		m(2, 8) = 2.0 * s[i] * (r[i] + s[i] - 1.0);
		m(2, 9) = (-2.0) * r[i] * (r[i] + s[i] - 1.0);
		m(2, 10) = 2.0 * r[i] * s[i];
		m(2, 11) = (-2.0) * s[i] * (r[i] + s[i] - 1.0);
		m(2, 12) = 2.0 * t[i] * (r[i] + s[i] - 1.0);
		m(2, 13) = (-2.0) * r[i] * t[i];
		m(2, 14) = (-2.0) * s[i] * t[i];
	}
}

void Prisma15::setBaseFunctions(Element &self)
{
	size_t GPCount = 9, nodeCount = 15;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(3, nodeCount));
	self.weighFactor = new std::vector<double>();

	std::vector< std::vector<double> > rst(3, std::vector<double>(GPCount));

	switch (GPCount) {
	case 9: {
		double v1 = 1.0 / 6.0;
		double v2 = 4.0 / 6.0;
		double v3 = sqrt(3.0 / 5.0);
		double v4 = 0.0;
		rst[0] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
		rst[1] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
		rst[2] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };
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
		m(0, 0) = -(1.0 - r - s) * (1.0 - t) * (2.0 * r + 2.0 * s + t) / 2.0;
		m(0, 1) = r * (1.0 - t) * (2.0 * r - t - 2.0) / 2.0;
		m(0, 2) = s * (1.0 - t) * (2.0 * s - t - 2.0) / 2.0;
		m(0, 3) = -(1.0 - r - s) * (1.0 + t) * (2.0 * r + 2.0 * s - t) / 2.0;
		m(0, 4) = r * (t + 1.0) * (2.0 * r + t - 2.0) / 2.0;
		m(0, 5) = s * (t + 1.0) * (2.0 * s + t - 2.0) / 2.0;
		m(0, 6) = 2.0 * r * (1.0 - r - s) * (1.0 - t);
		m(0, 7) = 2.0 * r * s * (1.0 - t);
		m(0, 8) = 2.0 * s * (1.0 - r - s) * (1.0 - t);
		m(0, 9) = 2.0 * r * (1.0 - r - s) * (1.0 + t);
		m(0, 10) = 2.0 * r * s * (1.0 + t);
		m(0, 11) = 2.0 * s * (1.0 - r - s) * (1.0 + t);
		m(0, 12) = (1.0 - r - s) * (1.0 - pow(t, 2.0));
		m(0, 13) = r * (1.0 - pow(t, 2.0));
		m(0, 14) = s * (1.0 - pow(t, 2.0));
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		///dN contains [dNr, dNs, dNt]
		MatrixDense &m = (*self.dN)[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// dNr - derivation of basis function
		m(0, 0) = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
		m(0, 1) = ((t - 1.0) * (t - 2.0 * r + 2.0)) / 2.0 - r * (t - 1.0);
		m(0, 2) = 0.0;
		m(0, 3) = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
		m(0, 4) = r * (t + 1.0) + ((t + 1.0) * (2.0 * r + t - 2.0)) / 2.0;
		m(0, 5) = 0.0;
		m(0, 6) = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * r * (t - 1.0);
		m(0, 7) = (-2.0) * s * (t - 1.0);
		m(0, 8) = 2.0 * s * (t - 1.0);
		m(0, 9) = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * r * (t + 1.0);
		m(0, 10) = 2.0 * s * (t + 1.0);
		m(0, 11) = -2.0 * s * (t + 1.0);
		m(0, 12) = pow(t, 2.0) - 1.0;
		m(0, 13) = 1.0 - pow(t, 2.0);
		m(0, 14) = 0.0;

		// dNs - derivation of basis function
		m(1, 0) = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
		m(1, 1) = 0.0;
		m(1, 2) = ((t - 1.0) * (t - 2.0 * s + 2.0)) / 2.0 - s * (t - 1.0);
		m(1, 3) = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
		m(1, 4) = 0.0;
		m(1, 5) = s * (t + 1.0) + ((t + 1.0) * (2.0 * s + t - 2.0)) / 2.0;
		m(1, 6) = 2.0 * r * (t - 1.0);
		m(1, 7) = (-2.0) * r * (t - 1.0);
		m(1, 8) = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * s * (t - 1.0);
		m(1, 9) = (-2.0) * r * (t + 1.0);
		m(1, 10) = 2.0 * r * (t + 1.0);
		m(1, 11) = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * s * (t + 1.0);
		m(1, 12) = pow(t, 2.0) - 1.0;
		m(1, 13) = 0.0;
		m(1, 14) = 1.0 - pow(t, 2.0);

		// dNt - derivation of basis function
		m(2, 0) = -((r + s - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0 - ((t - 1.0) * (r + s - 1.0)) / 2.0;
		m(2, 1) = (r * (t - 2.0 * r + 2.0)) / 2.0 + (r * (t - 1.0)) / 2.0;
		m(2, 2) = (s * (t - 2.0 * s + 2.0)) / 2.0 + (s * (t - 1.0)) / 2.0;
		m(2, 3) = ((r + s - 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 - ((t + 1.0) * (r + s - 1.0)) / 2.0;
		m(2, 4) = (r * (2.0 * r + t - 2.0)) / 2.0 + (r * (t + 1.0)) / 2.0;
		m(2, 5) = (s * (2.0 * s + t - 2.0)) / 2.0 + (s * (t + 1.0)) / 2.0;
		m(2, 6) = 2.0 * r * (r + s - 1.0);
		m(2, 7) = (-2.0) * r * s;
		m(2, 8) = 2.0 * s * (r + s - 1.0);
		m(2, 9) = (-2.0) * r * (r + s - 1.0);
		m(2, 10) = 2.0 * r * s;
		m(2, 11) = (-2.0) * s * (r + s - 1.0);
		m(2, 12) = 2.0 * t * (r + s - 1.0);
		m(2, 13) = (-2.0) * r * t;
		m(2, 14) = (-2.0) * s * t;
	}

	switch (GPCount) {
	case 9: {
		double v1 = 5.0 / 54.0;
		double v2 = 8.0 / 54.0;
		(*self.weighFactor) = { v1, v1, v1, v2, v2, v2, v1, v1, v1 };
		break;
	}
	default:
		exit(1);
	}

	BaseFunctions::created(self);
}


