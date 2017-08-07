
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS2D_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS2D_H_

#include "physics.h"
#include <cstring>

namespace espreso {

struct Physics2D: public virtual Physics {

	Physics2D() : Physics("", NULL, NULL) // skipped because Physics is inherited virtually
	{

	}

	double determinant2x2(double *values) const
	{
		double det = values[0] * values[3] - values[1] * values[2];
		if (det < 0) {
			printf("negative determinant\n");
			exit(0);
		}
		return det;
	}

	void inverse2x2(const double *m, double *inv, double det) const
	{
		double detJx = 1 / det;
		inv[0] =   detJx * m[3];
		inv[1] = - detJx * m[1];
		inv[2] = - detJx * m[2];
		inv[3] =   detJx * m[0];
	}

	// source dX, dY

	// target::
	// dX   0
	//  0  dY
	// dY  dX
	void distribute3x2(double *target, double *source, size_t rows, size_t columns) const
	{
		memcpy(target                               , source          , sizeof(double) * columns);
		memcpy(target + 2 * rows * columns + columns, source          , sizeof(double) * columns);

		memcpy(target + 1 * rows * columns + columns, source + columns, sizeof(double) * columns);
		memcpy(target + 2 * rows * columns          , source + columns, sizeof(double) * columns);
	}

	// source dX, dY

	// target::
	// dX   0
	//  0  dY
	//  0  0
	// dY  dX
	void distribute4x2(double *target, double *source, size_t rows, size_t columns) const
	{
		memcpy(target                               , source          , sizeof(double) * columns);
		memcpy(target + 3 * rows * columns + columns, source          , sizeof(double) * columns);

		memcpy(target + 1 * rows * columns + columns, source + columns, sizeof(double) * columns);
		memcpy(target + 3 * rows * columns          , source + columns, sizeof(double) * columns);
	}



	virtual void prepareHybridTotalFETIWithCorners();
	virtual void prepareHybridTotalFETIWithKernels();
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS2D_H_ */
