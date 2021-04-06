
#include "square.h"
#include "square4.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"
#include "basis/containers/point.h"

using namespace espreso;

void Square4::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, w;
	if (!Square::gpw(order, r, s, w)) {
		eslog::globalerror("ESPRESO internal error: cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(r.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(r.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		(*N)[i](0, 0) = 0.25 * (1 - r[i]) * (1 - s[i]);
		(*N)[i](0, 1) = 0.25 * (r[i] + 1) * (1 - s[i]);
		(*N)[i](0, 2) = 0.25 * (r[i] + 1) * (s[i] + 1);
		(*N)[i](0, 3) = 0.25 * (1 - r[i]) * (s[i] + 1);
	}

	for (size_t i = 0; i < r.size(); i++) {
		// dNs - derivation of basis function
		(*dN)[i](0, 0) = 0.25 * ( s[i] - 1);
		(*dN)[i](0, 1) = 0.25 * (-s[i] + 1);
		(*dN)[i](0, 2) = 0.25 * ( s[i] + 1);
		(*dN)[i](0, 3) = 0.25 * (-s[i] - 1);

		// dNt - derivation of basis function
		(*dN)[i](1, 0) = 0.25 * ( r[i] - 1);
		(*dN)[i](1, 1) = 0.25 * (-r[i] - 1);
		(*dN)[i](1, 2) = 0.25 * ( r[i] + 1);
		(*dN)[i](1, 3) = 0.25 * (-r[i] + 1);
	}

	BaseFunctions::created(*this);
}

void Square4::setBaseFunctions(Element &self)
{
	size_t GPCount = 4, nodeCount = 4;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(2, nodeCount));
	self.weighFactor = new std::vector<double>(GPCount, 1);

	{ // TODO: improve
		double CsQ_scale = 0.577350269189626;

		std::vector<double> s = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
		std::vector<double> t = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };

		for (unsigned int i = 0; i < GPCount; i++) {
			(*self.N)[i](0, 0) = 0.25 * (1 - s[i]) * (1 - t[i]);
			(*self.N)[i](0, 1) = 0.25 * (s[i] + 1) * (1 - t[i]);
			(*self.N)[i](0, 2) = 0.25 * (s[i] + 1) * (t[i] + 1);
			(*self.N)[i](0, 3) = 0.25 * (1 - s[i]) * (t[i] + 1);
		}

		for (unsigned int i = 0; i < GPCount; i++) {
			// dNs - derivation of basis function
			(*self.dN)[i](0, 0) = 0.25 * ( t[i] - 1);
			(*self.dN)[i](0, 1) = 0.25 * (-t[i] + 1);
			(*self.dN)[i](0, 2) = 0.25 * ( t[i] + 1);
			(*self.dN)[i](0, 3) = 0.25 * (-t[i] - 1);

			// dNt - derivation of basis function
			(*self.dN)[i](1, 0) = 0.25 * ( s[i] - 1);
			(*self.dN)[i](1, 1) = 0.25 * (-s[i] - 1);
			(*self.dN)[i](1, 2) = 0.25 * ( s[i] + 1);
			(*self.dN)[i](1, 3) = 0.25 * (-s[i] + 1);
		}
	}

	self.nN = new std::vector<MatrixDense>(nodeCount, MatrixDense(1, nodeCount));
	self.ndN = new std::vector<MatrixDense>(nodeCount, MatrixDense(2, nodeCount));

	{
		double CsQ_scale = 1;

		std::vector<double> s = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
		std::vector<double> t = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };

		for (unsigned int i = 0; i < nodeCount; i++) {
			(*self.nN)[i](0, 0) = 0.25 * (1 - s[i]) * (1 - t[i]);
			(*self.nN)[i](0, 1) = 0.25 * (s[i] + 1) * (1 - t[i]);
			(*self.nN)[i](0, 2) = 0.25 * (s[i] + 1) * (t[i] + 1);
			(*self.nN)[i](0, 3) = 0.25 * (1 - s[i]) * (t[i] + 1);
		}

		for (unsigned int i = 0; i < nodeCount; i++) {
			// node dNs - derivation of basis function
			(*self.ndN)[i](0, 0) = 0.25 * ( t[i] - 1);
			(*self.ndN)[i](0, 1) = 0.25 * (-t[i] + 1);
			(*self.ndN)[i](0, 2) = 0.25 * ( t[i] + 1);
			(*self.ndN)[i](0, 3) = 0.25 * (-t[i] - 1);

			// node dNt - derivation of basis function
			(*self.ndN)[i](1, 0) = 0.25 * ( s[i] - 1);
			(*self.ndN)[i](1, 1) = 0.25 * (-s[i] - 1);
			(*self.ndN)[i](1, 2) = 0.25 * ( s[i] + 1);
			(*self.ndN)[i](1, 3) = 0.25 * (-s[i] + 1);
		}
	}

	BaseFunctions::created(self);
}

void Square4::computeReferenceCoords(const MatrixDense & vertices, const MatrixDense & points, MatrixDense & result)
{
	MatrixDense tmp(4,2);
	tmp[0][0] = -vertices[0][0]+vertices[1][0]+vertices[2][0]-vertices[3][0];   tmp[0][1] = -vertices[0][1]+vertices[1][1]+vertices[2][1]-vertices[3][1]; // s
	tmp[1][0] = -vertices[0][0]-vertices[1][0]+vertices[2][0]+vertices[3][0];   tmp[1][1] = -vertices[0][1]-vertices[1][1]+vertices[2][1]+vertices[3][1]; // t
	tmp[2][0] =  vertices[0][0]-vertices[1][0]+vertices[2][0]-vertices[3][0];   tmp[2][1] =  vertices[0][1]-vertices[1][1]+vertices[2][1]-vertices[3][1]; // st
	tmp[3][0] =  vertices[0][0]+vertices[1][0]+vertices[2][0]+vertices[3][0];   tmp[3][1] =  vertices[0][1]+vertices[1][1]+vertices[2][1]+vertices[3][1]; // 1
	MatrixDense  F(1,2);
	MatrixDense dF(2,2);
	double bx, by, rx, ry, detdF, normF, detdF0 = tmp[0][0]*tmp[1][1]-tmp[0][1]*tmp[1][0];
	int cnt;
	result.resize(points.nrows,2);
	for (esint i = 0; i < points.nrows; ++i) {
		bx = 4.0*points[i][0]-tmp[3][0];
		by = 4.0*points[i][1]-tmp[3][1];
		rx = ( tmp[1][1]*bx-tmp[1][0]*by)/detdF0; // s
		ry = (-tmp[0][1]*bx+tmp[0][0]*by)/detdF0; // t
		cnt = 0;
		F[0][0]  = tmp[0][0]*rx+tmp[1][0]*ry+tmp[2][0]*rx*ry-bx;
		F[0][1]  = tmp[0][1]*rx+tmp[1][1]*ry+tmp[2][1]*rx*ry-by;
		dF[0][0] = tmp[0][0]+tmp[2][0]*ry;
		dF[0][1] = tmp[0][1]+tmp[2][1]*ry;
		dF[1][0] = tmp[1][0]+tmp[2][0]*rx;
		dF[1][1] = tmp[1][1]+tmp[2][1]*rx;
		detdF = tmp[0][0]*tmp[1][1]-tmp[0][1]*tmp[1][0];
		normF = std::sqrt(F[0][0]*F[0][0]+F[0][1]*F[0][1]);
		while (cnt > 0 && normF > 1e-5 && std::abs(detdF/detdF0) > 1e-5) {
			rx -= ( dF[1][1]*F[0][0]-dF[1][0]*F[0][1])/detdF;
			ry -= (-dF[1][0]*F[0][0]+dF[0][0]*F[0][1])/detdF;
			F[0][0]  = tmp[0][0]*rx+tmp[1][0]*ry+tmp[2][0]*rx*ry-bx;
			F[0][1]  = tmp[0][1]*rx+tmp[1][1]*ry+tmp[2][1]*rx*ry-by;
			dF[0][0] = tmp[0][0]+tmp[2][0]*ry;
			dF[0][1] = tmp[0][1]+tmp[2][1]*ry;
			dF[1][0] = tmp[1][0]+tmp[2][0]*rx;
			dF[1][1] = tmp[1][1]+tmp[2][1]*rx;
			detdF = tmp[0][0]*tmp[1][1]-tmp[0][1]*tmp[1][0];
			normF = std::sqrt(F[0][0]*F[0][0]+F[0][1]*F[0][1]);
			cnt--;
		}
		result[i][0] = rx;
		result[i][1] = ry;
	}
}
