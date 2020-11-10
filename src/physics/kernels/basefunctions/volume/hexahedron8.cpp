
#include "hexahedron.h"
#include "hexahedron8.h"


#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Hexahedron8::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, t, w;
	if (!Hexahedron::gpw(order, r, s, t, w)) {
		eslog::globalerror("ESPRESO internal error: cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(s.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(s.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		(*N)[i](0, 0) = 0.125 * (1 - r[i]) * (1 - s[i]) * (1 - t[i]);
		(*N)[i](0, 1) = 0.125 * (r[i] + 1) * (1 - s[i]) * (1 - t[i]);
		(*N)[i](0, 2) = 0.125 * (r[i] + 1) * (s[i] + 1) * (1 - t[i]);
		(*N)[i](0, 3) = 0.125 * (1 - r[i]) * (s[i] + 1) * (1 - t[i]);
		(*N)[i](0, 4) = 0.125 * (1 - r[i]) * (1 - s[i]) * (t[i] + 1);
		(*N)[i](0, 5) = 0.125 * (r[i] + 1) * (1 - s[i]) * (t[i] + 1);
		(*N)[i](0, 6) = 0.125 * (r[i] + 1) * (s[i] + 1) * (t[i] + 1);
		(*N)[i](0, 7) = 0.125 * (1 - r[i]) * (s[i] + 1) * (t[i] + 1);
	}

	for (size_t i = 0; i < r.size(); i++) {
		(*dN)[i](0, 0) = 0.125 * (-(1 - s[i]) * (1 - t[i]));
		(*dN)[i](0, 1) = 0.125 * ( (1 - s[i]) * (1 - t[i]));
		(*dN)[i](0, 2) = 0.125 * ( (1 + s[i]) * (1 - t[i]));
		(*dN)[i](0, 3) = 0.125 * (-(1 + s[i]) * (1 - t[i]));
		(*dN)[i](0, 4) = 0.125 * (-(1 - s[i]) * (1 + t[i]));
		(*dN)[i](0, 5) = 0.125 * ( (1 - s[i]) * (1 + t[i]));
		(*dN)[i](0, 6) = 0.125 * ( (1 + s[i]) * (1 + t[i]));
		(*dN)[i](0, 7) = 0.125 * (-(1 + s[i]) * (1 + t[i]));

		(*dN)[i](1, 0) = 0.125 * (-(1 - r[i]) * (1 - t[i]));
		(*dN)[i](1, 1) = 0.125 * (-(1 + r[i]) * (1 - t[i]));
		(*dN)[i](1, 2) = 0.125 * ( (1 + r[i]) * (1 - t[i]));
		(*dN)[i](1, 3) = 0.125 * ( (1 - r[i]) * (1 - t[i]));
		(*dN)[i](1, 4) = 0.125 * (-(1 - r[i]) * (1 + t[i]));
		(*dN)[i](1, 5) = 0.125 * (-(1 + r[i]) * (1 + t[i]));
		(*dN)[i](1, 6) = 0.125 * ( (1 + r[i]) * (1 + t[i]));
		(*dN)[i](1, 7) = 0.125 * ( (1 - r[i]) * (1 + t[i]));

		(*dN)[i](2, 0) = 0.125 * (-(1 - r[i]) * (1 - s[i]));
		(*dN)[i](2, 1) = 0.125 * (-(1 + r[i]) * (1 - s[i]));
		(*dN)[i](2, 2) = 0.125 * (-(1 + r[i]) * (1 + s[i]));
		(*dN)[i](2, 3) = 0.125 * (-(1 - r[i]) * (1 + s[i]));
		(*dN)[i](2, 4) = 0.125 * ( (1 - r[i]) * (1 - s[i]));
		(*dN)[i](2, 5) = 0.125 * ( (1 + r[i]) * (1 - s[i]));
		(*dN)[i](2, 6) = 0.125 * ( (1 + r[i]) * (1 + s[i]));
		(*dN)[i](2, 7) = 0.125 * ( (1 - r[i]) * (1 + s[i]));
	}
}

void Hexahedron8::setBaseFunctions(Element &self)
{
	size_t GPCount = 8, nodeCount = 8;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(3, nodeCount));
	self.weighFactor = new std::vector<double>(GPCount, 1);

	double CsQ_scale = 1 / std::sqrt(3);

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		(*self.N)[i](0, 0) = 0.125 * (1 - r) * (1 - s) * (1 - t);
		(*self.N)[i](0, 1) = 0.125 * (r + 1) * (1 - s) * (1 - t);
		(*self.N)[i](0, 2) = 0.125 * (r + 1) * (s + 1) * (1 - t);
		(*self.N)[i](0, 3) = 0.125 * (1 - r) * (s + 1) * (1 - t);
		(*self.N)[i](0, 4) = 0.125 * (1 - r) * (1 - s) * (t + 1);
		(*self.N)[i](0, 5) = 0.125 * (r + 1) * (1 - s) * (t + 1);
		(*self.N)[i](0, 6) = 0.125 * (r + 1) * (s + 1) * (t + 1);
		(*self.N)[i](0, 7) = 0.125 * (1 - r) * (s + 1) * (t + 1);
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// dNr - derivation of basis function
		(*self.dN)[i](0, 0) = 0.125 * (-(1 - s) * (1 - t));
		(*self.dN)[i](0, 1) = 0.125 * ( (1 - s) * (1 - t));
		(*self.dN)[i](0, 2) = 0.125 * ( (1 + s) * (1 - t));
		(*self.dN)[i](0, 3) = 0.125 * (-(1 + s) * (1 - t));
		(*self.dN)[i](0, 4) = 0.125 * (-(1 - s) * (1 + t));
		(*self.dN)[i](0, 5) = 0.125 * ( (1 - s) * (1 + t));
		(*self.dN)[i](0, 6) = 0.125 * ( (1 + s) * (1 + t));
		(*self.dN)[i](0, 7) = 0.125 * (-(1 + s) * (1 + t));

		// dNs - derivation of basis function
		(*self.dN)[i](1, 0) = 0.125 * (-(1 - r) * (1 - t));
		(*self.dN)[i](1, 1) = 0.125 * (-(1 + r) * (1 - t));
		(*self.dN)[i](1, 2) = 0.125 * ( (1 + r) * (1 - t));
		(*self.dN)[i](1, 3) = 0.125 * ( (1 - r) * (1 - t));
		(*self.dN)[i](1, 4) = 0.125 * (-(1 - r) * (1 + t));
		(*self.dN)[i](1, 5) = 0.125 * (-(1 + r) * (1 + t));
		(*self.dN)[i](1, 6) = 0.125 * ( (1 + r) * (1 + t));
		(*self.dN)[i](1, 7) = 0.125 * ( (1 - r) * (1 + t));

		// dNt - derivation of basis function
		(*self.dN)[i](2, 0) = 0.125 * (-(1 - r) * (1 - s));
		(*self.dN)[i](2, 1) = 0.125 * (-(1 + r) * (1 - s));
		(*self.dN)[i](2, 2) = 0.125 * (-(1 + r) * (1 + s));
		(*self.dN)[i](2, 3) = 0.125 * (-(1 - r) * (1 + s));
		(*self.dN)[i](2, 4) = 0.125 * ( (1 - r) * (1 - s));
		(*self.dN)[i](2, 5) = 0.125 * ( (1 + r) * (1 - s));
		(*self.dN)[i](2, 6) = 0.125 * ( (1 + r) * (1 + s));
		(*self.dN)[i](2, 7) = 0.125 * ( (1 - r) * (1 + s));
	}

	BaseFunctions::created(self);
}




