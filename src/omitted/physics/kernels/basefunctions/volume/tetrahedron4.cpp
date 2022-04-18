
#include "tetrahedron.h"
#include "tetrahedron4.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Tetrahedron4::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, t, w;
	if (!Tetrahedron::gpw(order, r, s, t, w)) {
		eslog::internalFailure("cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(s.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(s.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*N)[i];

		m(0, 0) = r[i];
		m(0, 1) = s[i];
		m(0, 2) = t[i];
		m(0, 3) = 1.0 - r[i] - s[i] - t[i];
	}

	for (size_t i = 0; i < r.size(); i++) {
		//  N = [ r, s, t,  1 - r - s - t ];
		MatrixDense &m = (*dN)[i];

		// dNr = [ 1, 0, 0, -1 ];
		m(0, 0) =  1.0;
		m(0, 1) =  0.0;
		m(0, 2) =  0.0;
		m(0, 3) = -1.0;

		// dNs = [ 0, 1, 0, -1 ];
		m(2, 0) =  0.0;
		m(2, 1) =  1.0;
		m(2, 2) =  0.0;
		m(2, 3) = -1.0;

		// dNs = [ 0, 0, 1, -1 ];
		m(1, 0) =  0.0;
		m(1, 1) =  0.0;
		m(1, 2) =  1.0;
		m(1, 3) = -1.0;
	}
}

void Tetrahedron4::setBaseFunctions(Element &self)
{
	size_t GPCount = 4, nodeCount = 4;

	std::vector<std::vector<double> > rst(3);

	switch (GPCount) {
	case 4: {
		rst[0] = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		rst[1] = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		rst[2] = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
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

		m(0, 0) = r;
		m(0, 1) = s;
		m(0, 2) = t;
		m(0, 3) = 1.0 - r - s - t;
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		//  N = [ r, s, t,  1 - r - s - t ];
		MatrixDense &m = (*self.dN)[i];

		// dNr = [ 1, 0, 0, -1 ];
		m(0, 0) =  1.0;
		m(0, 1) =  0.0;
		m(0, 2) =  0.0;
		m(0, 3) = -1.0;

		// dNs = [ 0, 1, 0, -1 ];
		m(2, 0) =  0.0;
		m(2, 1) =  1.0;
		m(2, 2) =  0.0;
		m(2, 3) = -1.0;

		// dNs = [ 0, 0, 1, -1 ];
		m(1, 0) =  0.0;
		m(1, 1) =  0.0;
		m(1, 2) =  1.0;
		m(1, 3) = -1.0;
	}

	switch (GPCount) {
	case 4: {
		self.weighFactor->resize(4, 1.0 / 24.0);
		break;
	}
	default:
		exit(1);
	}

	BaseFunctions::created(self);
}






