
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS3D_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS3D_H_

#include "physics.h"

namespace espreso {

struct Physics3D: public virtual Physics {

	Physics3D() : Physics("", NULL, NULL) // skipped because Physics is inherited virtually
	{

	}

	double determinant3x3(double *values) const
	{
		double det =
			+ values[0] * values[4] * values[8]
			+ values[1] * values[5] * values[6]
			+ values[2] * values[3] * values[7]
			- values[2] * values[4] * values[6]
			- values[1] * values[3] * values[8]
			- values[0] * values[5] * values[7];
		if (det < 0) {
			printf("nonpositive determinant\n");
			exit(0);
		}
		return det;
	}

	void inverse3x3(const double *m, double *inv, double det) const
	{
		double detJx = 1 / det;
		inv[0] = detJx * ( m[8] * m[4] - m[7] * m[5]);
		inv[1] = detJx * (-m[8] * m[1] + m[7] * m[2]);
		inv[2] = detJx * ( m[5] * m[1] - m[4] * m[2]);
		inv[3] = detJx * (-m[8] * m[3] + m[6] * m[5]);
		inv[4] = detJx * ( m[8] * m[0] - m[6] * m[2]);
		inv[5] = detJx * (-m[5] * m[0] + m[3] * m[2]);
		inv[6] = detJx * ( m[7] * m[3] - m[6] * m[4]);
		inv[7] = detJx * (-m[7] * m[0] + m[6] * m[1]);
		inv[8] = detJx * ( m[4] * m[0] - m[3] * m[1]);
	}

	void prepareHybridTotalFETIWithCorners();
	void prepareHybridTotalFETIWithKernels();
	void assembleB0FromCorners(const Step &step);
	void assembleB0FromKernels(const Step &step);
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS3D_H_ */
