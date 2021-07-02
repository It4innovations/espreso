
#ifndef SRC_ANALYSIS_MODE_NEWTONRAPHSON_H_
#define SRC_ANALYSIS_MODE_NEWTONRAPHSON_H_

#include "analysis/linearsystem/linearsystem.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct AX_NewtonRaphson {

	void init(AX_LinearSystem<double> *system);

private:
	AX_LinearSystem<double> *system;

	Matrix_Base<double> *K;
	Vector_Base<double> *U, *R, *f, *BC;
	Vector_Base<double> *lsSolution, *lsRHS, *lsResidual;
};

}

#endif /* SRC_ANALYSIS_MODE_NEWTONRAPHSON_H_ */
