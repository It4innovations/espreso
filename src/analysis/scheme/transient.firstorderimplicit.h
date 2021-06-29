
#ifndef SRC_ANALYSIS_SCHEME_TRANSIENT_FIRSTORDERIMPLICIT_H_
#define SRC_ANALYSIS_SCHEME_TRANSIENT_FIRSTORDERIMPLICIT_H_

#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct AX_TransientFirstOrderImplicit {

	void init(DirectSolver<double> *solver);
	void init(FETISolver<double> *solver);
	void init(MultigridSolver<double> *solver);

	Matrix_Base<double> *K, *M;
	Vector_Base<double> *f;

	Vector_Base<double> *U, *dU, *V, *X, *Y, *dTK, *dTM;
};

}

#endif /* SRC_ANALYSIS_SCHEME_TRANSIENT_FIRSTORDERIMPLICIT_H_ */
