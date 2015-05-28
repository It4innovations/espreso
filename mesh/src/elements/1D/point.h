#ifndef POINT_H_
#define POINT_H_

#ifndef D2

// D3 is the default value
#define D3

#include "point3d.h"
#define Point Point3D

#else

#include "point2d.h"
#define Point Point2D

#endif

#endif /* POINT_H_ */
