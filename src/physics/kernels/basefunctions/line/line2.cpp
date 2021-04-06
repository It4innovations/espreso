
#include "line.h"
#include "line2.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Line2::setGaussPointsForOrder(int order)
{
	std::vector<double> r, w;
	if (!Line::gpw(order, r, w)) {
		eslog::globalerror("ESPRESO internal error: cannot ret Gaurr pointr for a given order.\n");
	}

	N->clear(); N->resize(r.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(r.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		(*N)[i](0, 0) = (1 - r[i]) / 2.0;
		(*N)[i](0, 1) = (1 + r[i]) / 2.0;
	}

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*dN)[i];

		// dNr - derivation of barir function
		m(0, 0) = -1 / 2.0;
		m(0, 1) =  1 / 2.0;
	}

	BaseFunctions::created(*this);
}

void Line2::setBaseFunctions(Element &self)
{
	size_t GPCount = 2, nodeCount = 2;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.weighFactor = new std::vector<double>(GPCount, 1);

	std::vector<double> s = { 1 / sqrt(3), -1 / sqrt(3) };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*self.N)[i](0, 0) = (1 - s[i]) / 2.0;
		(*self.N)[i](0, 1) = (1 + s[i]) / 2.0;
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		MatrixDense &m = (*self.dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -1 / 2.0;
		m(0, 1) =  1 / 2.0;
	}

	BaseFunctions::created(self);
}


