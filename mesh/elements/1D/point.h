#ifndef POINT_H_
#define POINT_H_

#define ESPRESO_POINT_DIMENSION 3

#if ESPRESO_POINT_DIMENSION == 2
	#include "point2d.h"

	namespace mesh {
		typedef mesh::Point2D Point;
	}
#elif ESPRESO_POINT_DIMENSION == 3
	#include "point3d.h"

	namespace mesh {
		typedef mesh::Point3D Point;
	}
#else
	#error "Incorrect user-supplied value for ESPRESO_POINT_DIMENSION"
#endif

#endif /* POINT_H_ */
