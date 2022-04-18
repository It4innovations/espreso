
#include "basefunctions.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "mesh/mesh.h"
#include "math/matrix.dense.h"
#include "math/math.h"

#include "point/point1.h"
#include "line/line2.h"
#include "line/line3.h"
#include "plane/square4.h"
#include "plane/square8.h"
#include "plane/triangle3.h"
#include "plane/triangle6.h"
#include "volume/hexahedron20.h"
#include "volume/hexahedron8.h"
#include "volume/prisma15.h"
#include "volume/prisma6.h"
#include "volume/pyramid13.h"
#include "volume/pyramid5.h"
#include "volume/tetrahedron10.h"
#include "volume/tetrahedron4.h"

using namespace espreso;

void BaseFunctions::setBaseFunctions()
{
	Point1       ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::POINT1   )]);
	Line2        ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::LINE2    )]);
	Triangle3    ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE3)]);
	Square4      ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::SQUARE4  )]);
	Tetrahedron4 ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::TETRA4   )]);
	Pyramid5     ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::PYRAMID5 )]);
	Prisma6      ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::PRISMA6  )]);
	Hexahedron8  ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::HEXA8    )]);

	Line3        ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::LINE3    )]);
	Triangle6    ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE6)]);
	Square8      ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::SQUARE8  )]);
	Tetrahedron10::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::TETRA10  )]);
	Pyramid13    ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::PYRAMID13)]);
	Prisma15     ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::PRISMA15 )]);
	Hexahedron20 ::setBaseFunctions(Mesh::edata[static_cast<int>(Element::CODE::HEXA20   )]);
}

BaseFunctions::~BaseFunctions()
{
	if (N != NULL) { delete N; }
	if (NN != NULL) { delete NN; }
	if (NNN != NULL) { delete NNN; }
	if (dN != NULL) { delete dN; }
	if (weighFactor != NULL) { delete weighFactor; }

	if (nN != NULL) { delete nN; }
	if (ndN != NULL) { delete ndN; }
}

void BaseFunctions::created(Element &e)
{
	if (e.N) {
		if (e.NN) { delete e.NN; }
		if (e.NNN) { delete e.NNN; }
		e.NN  = new std::vector<MatrixDense>(e.N->size(), MatrixDense(2 * e.N->front().nrows, 2 * e.N->front().ncols));
		e.NNN = new std::vector<MatrixDense>(e.N->size(), MatrixDense(3 * e.N->front().nrows, 3 * e.N->front().ncols));
		for (size_t gp = 0; gp < e.N->size(); gp++) {
			e.NN->at(gp).fill(0);
			e.NNN->at(gp).fill(0);
			for (esint r = 0; r < e.N->front().nrows; r++) {
				for (esint c = 0; c < e.N->front().ncols; c++) {
					e.NN->at(gp)( 0 * e.N->front().nrows + r, 0 * e.N->front().ncols + c) = e.N->at(gp)(r, c);
					e.NN->at(gp)( 1 * e.N->front().nrows + r, 1 * e.N->front().ncols + c) = e.N->at(gp)(r, c);
					e.NNN->at(gp)(0 * e.N->front().nrows + r, 0 * e.N->front().ncols + c) = e.N->at(gp)(r, c);
					e.NNN->at(gp)(1 * e.N->front().nrows + r, 1 * e.N->front().ncols + c) = e.N->at(gp)(r, c);
					e.NNN->at(gp)(2 * e.N->front().nrows + r, 2 * e.N->front().ncols + c) = e.N->at(gp)(r, c);
				}
			}
		}
	}
}

void BaseFunctions::recomputeDetJ(Element *e, MatrixDense& coords, MatrixDense& resdetJ, MatrixDense* points)
{
	MatrixDense J1(1, 1), J2(2, 3), J3(3, 3);
	double detJ = 0;
	if (points == NULL) {
		resdetJ.resize(1, e->weighFactor->size());
		for (size_t gp = 0; gp < e->weighFactor->size(); gp++) {
			switch (e->type) {
			case Element::TYPE::POINT:
				detJ = 1.0;
				break;
			case Element::TYPE::LINE:
				J1.multiply((*e->dN)[gp], coords);
				detJ = J1.vals[0];
				break;
			case Element::TYPE::PLANE:
				J2.multiply((*e->dN)[gp], coords);
				detJ = J2[0][0]*J2[1][1] - J2[0][1]*J2[1][0];
				break;
			case Element::TYPE::VOLUME:
				J3.multiply((*e->dN)[gp], coords);
				detJ = MATH::determinant3x3(J3.vals);
				break;
			}
			resdetJ[0][gp] = detJ;
		}
	} else {
		resdetJ.resize(1,points->nrows);
		for (esint p = 0; p < points->nrows; p++) {
			switch (e->type) {
			case Element::TYPE::POINT:
				detJ = 1.0;
				break;
			case Element::TYPE::LINE:
				J1.multiply((*e->dN)[p], coords);
				detJ = J1.vals[0];
				break;
			case Element::TYPE::PLANE:
				J2.multiply((*e->dN)[p], coords);
				detJ = J2[0][0]*J2[1][1] - J2[0][1]*J2[1][0];
				break;
			case Element::TYPE::VOLUME:
				J3.multiply((*e->dN)[p], coords);
				detJ = MATH::determinant3x3(J3.vals);
				break;
			}
			resdetJ[0][p] = detJ;
		}
	}
}

void BaseFunctions::recomputeDetJN(Element *e, MatrixDense& coords, MatrixDense& resdetJ, MatrixDense& resN, MatrixDense& refPoints)
{
	MatrixDense J1(1, 1), J2(2, 2), J3(3, 3), t_dN;
	double detJ = 0;
	resdetJ.resize(1,refPoints.nrows);
	resN.resize(e->nodes, refPoints.nrows);
	t_dN.resize(2, e->nodes);
	for (esint p = 0; p < refPoints.nrows; p++) {
		switch (e->code) {
		case Element::CODE::TRIANGLE3:
			resN[0][p] =  1 - refPoints[p][0] - refPoints[p][1];
			resN[1][p] = refPoints[p][0];
			resN[2][p] = refPoints[p][1];
			t_dN[0][0] = -1;   t_dN[1][0] = -1;
			t_dN[0][1] =  1;   t_dN[1][1] =  0;
			t_dN[0][2] =  0;   t_dN[1][2] =  1;
			break;
		case Element::CODE::SQUARE4:
			resN[0][p] = 0.25 * (1 - refPoints[p][0]) * (1 - refPoints[p][1]);
			resN[1][p] = 0.25 * (1 + refPoints[p][0]) * (1 - refPoints[p][1]);
			resN[2][p] = 0.25 * (1 + refPoints[p][0]) * (1 + refPoints[p][1]);
			resN[3][p] = 0.25 * (1 - refPoints[p][0]) * (1 + refPoints[p][1]);
			t_dN[0][0] = 0.25 * ( refPoints[p][1] - 1);   t_dN[1][0] = 0.25 * ( refPoints[p][0] - 1);
			t_dN[0][1] = 0.25 * (-refPoints[p][1] + 1);   t_dN[1][1] = 0.25 * (-refPoints[p][0] - 1);
			t_dN[0][2] = 0.25 * ( refPoints[p][1] + 1);   t_dN[1][2] = 0.25 * ( refPoints[p][0] + 1);
			t_dN[0][3] = 0.25 * (-refPoints[p][1] - 1);   t_dN[1][3] = 0.25 * (-refPoints[p][0] + 1);
			break;
		default:
			break;
		}
		switch (e->type) {
		case Element::TYPE::POINT:
			detJ = 1.0;
			break;
		case Element::TYPE::LINE:
			J1.multiply(t_dN, coords);
			detJ = J1.vals[0];
			break;
		case Element::TYPE::PLANE:
			J2.multiply(t_dN, coords);
			detJ = J2[0][0]*J2[1][1] - J2[0][1]*J2[1][0];
			break;
		case Element::TYPE::VOLUME:
			J3.multiply(t_dN, coords);
			detJ = MATH::determinant3x3(J3.vals);
			break;
		}
		resdetJ[0][p] = detJ;
	}
}
