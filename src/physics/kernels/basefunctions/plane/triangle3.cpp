
#include "triangle.h"
#include "triangle3.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "math/matrix.dense.h"

using namespace espreso;

void Triangle3::setGaussPointsForOrder(int order)
{
	std::vector<double> r, s, w;
	if (!Triangle::gpw(order, r, s, w)) {
		eslog::globalerror("ESPRESO internal error: cannot set Gauss points for a given order.\n");
	}

	N->clear(); N->resize(r.size(), MatrixDense(1, nodes));
	dN->clear(); dN->resize(r.size(), MatrixDense(1, nodes));
	weighFactor->assign(w.begin(), w.end());

	for (size_t i = 0; i < r.size(); i++) {
		(*N)[i](0, 0) = 1 - r[i] - s[i];
		(*N)[i](0, 1) = r[i];
		(*N)[i](0, 2) = s[i];
	}

	for (size_t i = 0; i < r.size(); i++) {
		MatrixDense &m = (*dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -1;
		m(0, 1) =  1;
		m(0, 2) =  0;

		// dNt - derivation of basis function
		m(1, 0) = -1;
		m(1, 1) =  0;
		m(1, 2) =  1;
	}
}

void Triangle3::setBaseFunctions(Element &self)
{
	size_t GPCount = 6, nodeCount = 3;
	self.weighFactor = new std::vector<double>({    0.111690794839005, 0.111690794839005, 0.111690794839005, 0.054975871827661, 0.054975871827661, 0.054975871827661});
	std::vector<double> s = {                  0.445948490915965, 0.445948490915965, 0.108103018168070, 0.091576213509771, 0.091576213509771, 0.816847572980459};
	std::vector<double> t = {                  0.445948490915965, 0.108103018168070, 0.445948490915965, 0.091576213509771, 0.816847572980459, 0.091576213509771};

//	size_t GPCount = 1, nodeCount = 3;

	self.N = new std::vector<MatrixDense>(GPCount, MatrixDense(1, nodeCount));
	self.dN = new std::vector<MatrixDense>(GPCount, MatrixDense(2, nodeCount));
//	weighFactor = new std::vector<double>({ 1 / 2.0 });
//
//	std::vector<double> s = { 1.0 / 3 };
//	std::vector<double> t = { 1.0 / 3 };

	for (unsigned int i = 0; i < GPCount; i++) {
		(*self.N)[i](0, 0) = 1 - s[i] - t[i];
		(*self.N)[i](0, 1) = s[i];
		(*self.N)[i](0, 2) = t[i];
	}

	for (unsigned int i = 0; i < GPCount; i++) {
		MatrixDense &m = (*self.dN)[i];

		// dNs - derivation of basis function
		m(0, 0) = -1;
		m(0, 1) =  1;
		m(0, 2) =  0;

		// dNt - derivation of basis function
		m(1, 0) = -1;
		m(1, 1) =  0;
		m(1, 2) =  1;
	}

	BaseFunctions::created(self);
}

void Triangle3::computeReferenceCoords(const MatrixDense & vertices, const MatrixDense & points, MatrixDense & result)
{
	// assume non-zero denominator
	result.resize(points.nrows, 2);

	double ux = vertices[1][0] - vertices[0][0], uy = vertices[1][1] - vertices[0][1];
	double vx = vertices[2][0] - vertices[0][0], vy = vertices[2][1] - vertices[0][1];
	double uv = ux * vx + uy * vy, uu = ux * ux + uy * uy, vv = vx * vx + vy * vy;
	double denominator = uv * uv - uu * vv;
	for (esint r = 0; r < points.nrows; ++r) {
		double wx = points[r][0] - vertices[0][0], wy = points[r][1] - vertices[0][1];
		double wu = wx * ux + wy * uy, wv = wx * vx + wy * vy;
		result[r][0] = (uv * wv - vv * wu) / denominator;
		result[r][1] = (uv * wu - uu * wv) / denominator;
	}
}




