// Dalibor Lukas, March 2009


#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "bem.func.h"

extern int LineQuadratureOrders[];
extern int LineQuadratureSizes[];
extern double *LineQuadraturePoints[];
extern double *LineQuadratureWeights[];

extern int TriangleQuadratureOrders[];
extern int TriangleQuadratureSizes[];
extern Point2D *TriangleQuadraturePoints[];
extern double *TriangleQuadratureWeights[];

/*
extern int *RectangleQuadratureOrders;
extern int *RectangleQuadratureSizes;
extern double **RectangleQuadraturePoints;
extern double **RectangleQuadratureWeights;
*/

extern int TetrahedronQuadratureOrders[];
extern int TetrahedronQuadratureSizes[];
extern Point3D *TetrahedronQuadraturePoints[];
extern double *TetrahedronQuadratureWeights[];

#endif
