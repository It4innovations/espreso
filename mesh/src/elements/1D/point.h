#ifndef POINT_H_
#define POINT_H_

#if ESPRESO_POINT_DIMENSION == 2
	#define D2
	#include "point2d.h"

	namespace mesh {
		typedef mesh::Point2D Point;
	}
#else
	#define D3
	#include "point3d.h"

	namespace mesh {
		typedef mesh::Point3D Point;
	}
#endif

#endif /* POINT_H_ */
