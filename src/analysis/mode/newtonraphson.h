
#ifndef SRC_ANALYSIS_MODE_NEWTONRAPHSON_H_
#define SRC_ANALYSIS_MODE_NEWTONRAPHSON_H_

#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct AX_NewtonRaphson {

	void init(DirectSolver<double> *solver);
	void init(FETISolver<double> *solver);
	void init(MultigridSolver<double> *solver);

private:
	LinearSolver *solver;

	Matrix_Base<double> *K;
	Vector_Base<double> *U, *R, *f, *BC;
	Vector_Base<double> *lsSolution, *lsRHS, *lsResidual;
};

}

#endif /* SRC_ANALYSIS_MODE_NEWTONRAPHSON_H_ */
