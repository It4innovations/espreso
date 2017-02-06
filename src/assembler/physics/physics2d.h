
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS2D_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS2D_H_

#include "physics.h"

namespace espreso {

struct Physics2D: public Physics {

	Physics2D(Mesh *mesh, Instance *instance);

	double determinant2x2(double *values) const
	{
		return values[0] * values[3] - values[1] * values[2];
	}

	void inverse2x2(const double *m, double *inv, double det) const
	{
		double detJx = 1 / det;
		inv[0] =   detJx * m[3];
		inv[1] = - detJx * m[1];
		inv[2] = - detJx * m[2];
		inv[3] =   detJx * m[0];
	}

	void prepareHybridTotalFETIWithCorners();
	void prepareHybridTotalFETIWithKernels();
	void assembleB0FromCorners(const Step &step);
	void assembleB0FromKernels(const Step &step);
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS2D_H_ */
